from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar import cli_main
from drakkar import cli_parser
from drakkar import database_latest
from drakkar.database_latest import (
    LATEST_SOURCES,
    STATUS_AHEAD,
    STATUS_CURRENT,
    STATUS_OUTDATED,
    STATUS_UNCONFIGURED,
    STATUS_UNKNOWN,
    _latest_from_index,
    _latest_from_mtime,
    _latest_from_probe,
    _version_key,
    installed_version,
    normalize_latest_names,
    resolve_latest,
    run_database_latest,
)


PFAM_INDEX = """
<a href="Pfam36.0/">Pfam36.0/</a>
<a href="Pfam37.4/">Pfam37.4/</a>
<a href="Pfam38.2/">Pfam38.2/</a>
<a href="Pfam_SARS-CoV-2_2.0/">Pfam_SARS-CoV-2_2.0/</a>
"""

AMR_INDEX = """
<a href="2025-07-16.1/">2025-07-16.1/</a>
<a href="2026-08-07.1/">2026-08-07.1/</a>
<a href="latest/">latest/</a>
"""


class VersionKeyTests(unittest.TestCase):
    def test_pfam_minor_release_sorts_above_earlier_major(self) -> None:
        self.assertGreater(_version_key("Pfam38.2"), _version_key("Pfam37.4"))

    def test_kofam_archive_dates_sort_chronologically(self) -> None:
        self.assertGreater(_version_key("2026-07-02"), _version_key("2026-02-01"))

    def test_double_digit_dbcan_release_sorts_above_single_digit(self) -> None:
        self.assertGreater(_version_key("V15"), _version_key("V9"))

    def test_keys_of_differently_shaped_versions_compare_without_raising(self) -> None:
        self.assertNotEqual(_version_key("Pfam38.2"), _version_key("2026-07-02"))


class IndexDiscoveryTests(unittest.TestCase):
    def test_index_returns_newest_matching_release(self) -> None:
        with patch.object(database_latest, "_fetch_text", return_value=PFAM_INDEX):
            latest, _ = _latest_from_index(LATEST_SOURCES["pfam"], timeout=1)
        self.assertEqual(latest, "Pfam38.2")

    def test_index_ignores_releases_that_do_not_match_the_pattern(self) -> None:
        with patch.object(database_latest, "_fetch_text", return_value=AMR_INDEX):
            latest, _ = _latest_from_index(LATEST_SOURCES["amr"], timeout=1)
        # The NCBI directory carries a "latest/" alias that is not a release name.
        self.assertEqual(latest, "2026-08-07.1")

    def test_index_without_any_matching_release_raises(self) -> None:
        with patch.object(database_latest, "_fetch_text", return_value="<a href='none/'>none</a>"):
            with self.assertRaises(ValueError):
                _latest_from_index(LATEST_SOURCES["pfam"], timeout=1)


class ProbeDiscoveryTests(unittest.TestCase):
    def test_probe_walks_forward_until_a_release_is_missing(self) -> None:
        available = {"V14", "V15", "V16"}
        with patch.object(
            database_latest,
            "_url_exists",
            side_effect=lambda url, timeout: any(f"Databases/{v}/" in url for v in available),
        ):
            latest, _ = _latest_from_probe(LATEST_SOURCES["cazy"], "V14", timeout=1)
        self.assertEqual(latest, "V16")

    def test_probe_reports_the_installed_release_when_nothing_newer_exists(self) -> None:
        with patch.object(
            database_latest,
            "_url_exists",
            side_effect=lambda url, timeout: "Databases/V14/" in url,
        ):
            latest, _ = _latest_from_probe(LATEST_SOURCES["cazy"], "V14", timeout=1)
        self.assertEqual(latest, "V14")

    def test_probe_falls_back_to_the_floor_when_no_version_is_installed(self) -> None:
        seen = []

        def exists(url, timeout):
            seen.append(url)
            return "Databases/V14/" in url

        with patch.object(database_latest, "_url_exists", side_effect=exists):
            latest, _ = _latest_from_probe(LATEST_SOURCES["cazy"], None, timeout=1)
        self.assertEqual(latest, "V14")
        self.assertTrue(seen[0].endswith("dbCAN-HMMdb-V14.txt"))


class MtimeDiscoveryTests(unittest.TestCase):
    def test_last_modified_header_becomes_the_source_date(self) -> None:
        class Response:
            headers = {"Last-Modified": "Fri, 21 Aug 2026 11:30:21 GMT"}

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

        with patch.object(database_latest, "_open", return_value=Response()):
            latest, detail = _latest_from_mtime(LATEST_SOURCES["vfdb"], timeout=1)
        self.assertEqual(latest, "2026-08-21")
        self.assertIn("rolling", detail)

    def test_missing_last_modified_header_raises(self) -> None:
        class Response:
            headers: dict[str, str] = {}

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

        with patch.object(database_latest, "_open", return_value=Response()):
            with self.assertRaises(ValueError):
                _latest_from_mtime(LATEST_SOURCES["vfdb"], timeout=1)


class InstalledVersionTests(unittest.TestCase):
    def test_managed_database_version_comes_from_the_release_directory(self) -> None:
        self.assertEqual(
            installed_version("kegg", LATEST_SOURCES["kegg"], "/db/kofams/2026-02-01"),
            "2026-02-01",
        )

    def test_managed_database_version_prefers_the_recorded_version(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            release_dir = Path(directory) / "renamed_by_hand"
            release_dir.mkdir()
            (release_dir / "database_versions.yaml").write_text(
                "requested_version: 2026-02-01\n", encoding="utf-8"
            )
            self.assertEqual(
                installed_version("kegg", LATEST_SOURCES["kegg"], str(release_dir)),
                "2026-02-01",
            )

    def test_gtdb_version_is_read_out_of_the_release_directory_name(self) -> None:
        self.assertEqual(
            installed_version("gtdb", LATEST_SOURCES["gtdb"], "/db/gtdbtk_db/20260426/release232"),
            "232",
        )

    def test_version_is_none_when_the_path_does_not_carry_one(self) -> None:
        self.assertIsNone(
            installed_version("gtdb", LATEST_SOURCES["gtdb"], "/db/gtdbtk_db/20241001")
        )

    def test_version_is_none_when_nothing_is_configured(self) -> None:
        self.assertIsNone(installed_version("kegg", LATEST_SOURCES["kegg"], None))


class ResolveLatestTests(unittest.TestCase):
    def _resolve(self, name, latest, config):
        with patch.object(database_latest, "_latest_from_index", return_value=(latest, "")):
            return resolve_latest(name, config=config, timeout=1)

    def test_configured_release_matching_the_source_is_up_to_date(self) -> None:
        result = self._resolve("pfam", "Pfam37.4", {"PFAM_DB": "/db/pfam/Pfam37.4"})
        self.assertEqual(result.status, STATUS_CURRENT)

    def test_older_configured_release_is_outdated(self) -> None:
        result = self._resolve("pfam", "Pfam38.2", {"PFAM_DB": "/db/pfam/Pfam37.4"})
        self.assertEqual(result.status, STATUS_OUTDATED)
        self.assertEqual(result.installed, "Pfam37.4")
        self.assertEqual(result.latest, "Pfam38.2")

    def test_release_newer_than_the_source_listing_is_flagged_as_ahead(self) -> None:
        result = self._resolve("pfam", "Pfam37.4", {"PFAM_DB": "/db/pfam/Pfam38.2"})
        self.assertEqual(result.status, STATUS_AHEAD)

    def test_empty_config_key_is_reported_as_unconfigured(self) -> None:
        result = self._resolve("pfam", "Pfam38.2", {"PFAM_DB": "  "})
        self.assertEqual(result.status, STATUS_UNCONFIGURED)

    def test_unreachable_source_is_reported_as_unknown_without_raising(self) -> None:
        with patch.object(
            database_latest,
            "_latest_from_index",
            side_effect=urllib.error.URLError("host is down"),
        ):
            result = resolve_latest("pfam", config={"PFAM_DB": "/db/pfam/Pfam37.4"}, timeout=1)
        self.assertEqual(result.status, STATUS_UNKNOWN)
        self.assertIsNone(result.latest)
        self.assertIn("could not reach the source", result.detail)

    def test_outdated_managed_release_carries_a_reinstall_command(self) -> None:
        result = self._resolve("pfam", "Pfam38.2", {"PFAM_DB": "/db/pfam/Pfam37.4"})
        self.assertEqual(
            result.install_command,
            "drakkar database pfam --directory /db/pfam --version Pfam38.2",
        )

    def test_rolling_database_is_reinstalled_without_a_version(self) -> None:
        with patch.object(database_latest, "_latest_from_mtime", return_value=("2026-08-21", "")):
            result = resolve_latest("vfdb", config={"VFDB_DB": "/db/vfdb/2026-04-24"}, timeout=1)
        self.assertEqual(result.status, STATUS_OUTDATED)
        self.assertEqual(result.install_command, "drakkar database vfdb --directory /db/vfdb")

    def test_unconfigured_database_falls_back_to_the_base_directory(self) -> None:
        with patch.object(database_latest, "_latest_from_index", return_value=("Pfam38.2", "")):
            result = resolve_latest(
                "pfam", config={"PFAM_DB": "", "DATABASES_DIR": "/db"}, timeout=1
            )
        self.assertEqual(result.status, STATUS_UNCONFIGURED)
        self.assertEqual(
            result.install_command,
            "drakkar database pfam --directory /db/pfam --version Pfam38.2",
        )

    def test_unconfigured_database_without_a_base_directory_uses_a_placeholder(self) -> None:
        with patch.object(database_latest, "_latest_from_index", return_value=("Pfam38.2", "")):
            result = resolve_latest("pfam", config={"PFAM_DB": ""}, timeout=1)
        self.assertIn("<directory>", result.install_command)

    def test_unmanaged_database_reports_an_install_hint_instead_of_a_command(self) -> None:
        with patch.object(database_latest, "_latest_from_text", return_value=("232", "")):
            result = resolve_latest("gtdb", config={"GTDB_DB": "/db/gtdb/release226"}, timeout=1)
        self.assertEqual(result.status, STATUS_OUTDATED)
        self.assertIn("gtdbtk download-db.sh", result.install_command)


class NameNormalizationTests(unittest.TestCase):
    def test_no_names_selects_every_source(self) -> None:
        selected, unknown = normalize_latest_names([])
        self.assertEqual(selected, list(LATEST_SOURCES))
        self.assertEqual(unknown, [])

    def test_managed_alias_resolves_to_its_database(self) -> None:
        selected, unknown = normalize_latest_names(["kofams"])
        self.assertEqual(selected, ["kegg"])
        self.assertEqual(unknown, [])

    def test_repeated_names_are_selected_once(self) -> None:
        selected, _ = normalize_latest_names(["kegg", "kofams"])
        self.assertEqual(selected, ["kegg"])

    def test_unrecognised_name_is_reported(self) -> None:
        selected, unknown = normalize_latest_names(["kegg", "nonsense"])
        self.assertEqual(selected, ["kegg"])
        self.assertEqual(unknown, ["nonsense"])


class RunDatabaseLatestTests(unittest.TestCase):
    def test_unknown_database_name_fails_without_querying_any_source(self) -> None:
        with patch.object(database_latest, "resolve_all") as resolve_all:
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_latest(["nonsense"])
        self.assertEqual(exit_code, 1)
        resolve_all.assert_not_called()
        self.assertIn("Unknown database", output.getvalue())

    def test_outdated_database_is_reported_but_does_not_fail_the_command(self) -> None:
        with patch.object(database_latest, "_latest_from_index", return_value=("Pfam38.2", "")):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_latest(["pfam"], config={"PFAM_DB": "/db/pfam/Pfam37.4"})
        rendered = output.getvalue()
        self.assertEqual(exit_code, 0)
        self.assertIn("Pfam37.4", rendered)
        self.assertIn("Pfam38.2", rendered)
        self.assertIn("outdated", rendered)

    def test_unreachable_source_does_not_fail_the_command(self) -> None:
        with patch.object(
            database_latest,
            "_latest_from_index",
            side_effect=OSError("timed out"),
        ):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_latest(["pfam"], config={"PFAM_DB": "/db/pfam/Pfam37.4"})
        self.assertEqual(exit_code, 0)
        self.assertIn("unknown", output.getvalue())


class DatabaseLatestDispatchTests(unittest.TestCase):
    def test_latest_is_parsed_as_a_database_subcommand(self) -> None:
        args = cli_parser.build_parser().parse_args(["database", "latest", "kegg"])
        self.assertEqual(args.database_name, "latest")
        self.assertEqual(args.databases, ["kegg"])

    def test_latest_dispatches_without_launching_snakemake_or_checking_screen(self) -> None:
        argv = ["drakkar", "database", "latest", "kegg", "--timeout", "5"]

        with patch.object(cli_module, "run_database_latest", return_value=0) as checked, \
                patch.object(cli_module, "run_snakemake_database") as launched, \
                patch.object(cli_main, "check_screen_session") as screen, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 0)
        launched.assert_not_called()
        screen.assert_not_called()
        checked.assert_called_once_with(["kegg"], timeout=5)

    def test_installing_a_database_still_reaches_snakemake(self) -> None:
        args = cli_parser.build_parser().parse_args(
            ["database", "kegg", "--directory", "/db/kofams", "--version", "2026-02-01"]
        )
        self.assertEqual(args.database_name, "kegg")
        self.assertFalse(hasattr(args, "databases"))


if __name__ == "__main__":
    unittest.main()
