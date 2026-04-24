from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.cli import (
    replace_config_value,
    validate_database_version,
    validate_download_runtime,
    validate_managed_database_version,
)
from drakkar.database_registry import (
    database_artifact_path,
    database_release_from_config,
    database_source_version_label,
    database_sources,
    database_release_dir,
    database_target_path,
    normalize_managed_database_name,
)


class DatabaseCommandTests(unittest.TestCase):
    def test_normalize_managed_database_name_accepts_database_name(self) -> None:
        self.assertEqual(normalize_managed_database_name("amr"), "amr")

    def test_normalize_managed_database_name_accepts_kofams_alias(self) -> None:
        self.assertEqual(normalize_managed_database_name("kofams"), "kegg")

    def test_validate_database_version_rejects_paths(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNone(validate_database_version("2026/04/21"))

    def test_validate_download_runtime_accepts_positive_integer(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(validate_download_runtime(120), 120)

    def test_validate_download_runtime_rejects_non_positive_integer(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNone(validate_download_runtime(0))

    def test_validate_managed_database_version_accepts_kegg_archive_date(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(validate_managed_database_version("kegg", "2026-02-01"), "2026-02-01")

    def test_validate_managed_database_version_rejects_invalid_kegg_archive_date(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNone(validate_managed_database_version("kegg", "20260201"))

    def test_validate_managed_database_version_defaults_vfdb_to_download_date(self) -> None:
        with patch.object(cli_module, "default_database_version", return_value="2026-04-24"):
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(validate_managed_database_version("vfdb", None), "2026-04-24")

    def test_database_release_dir_joins_base_and_version(self) -> None:
        release_dir = database_release_dir("amr", "/tmp/amr", "20260421")
        self.assertEqual(release_dir, Path("/tmp/amr/20260421"))

    def test_database_target_path_uses_database_specific_basename(self) -> None:
        target_path = database_target_path("kegg", "/tmp/kegg", "20260421")
        self.assertEqual(target_path, Path("/tmp/kegg/20260421/kofams"))

    def test_database_release_from_config_accepts_release_directory(self) -> None:
        self.assertEqual(
            database_release_from_config("pfam", "/tmp/pfam/Pfam37.4"),
            Path("/tmp/pfam/Pfam37.4"),
        )

    def test_database_artifact_path_accepts_release_directory(self) -> None:
        self.assertEqual(
            database_artifact_path("pfam", "/tmp/pfam/Pfam37.4", "_ec.tsv"),
            Path("/tmp/pfam/Pfam37.4/pfam_ec.tsv"),
        )

    def test_database_artifact_path_accepts_legacy_primary_path(self) -> None:
        self.assertEqual(
            database_artifact_path("amr", "/tmp/amr/2025-07-16.1/amr", ".tsv"),
            Path("/tmp/amr/2025-07-16.1/amr.tsv"),
        )

    def test_kegg_database_sources_use_requested_archive_version(self) -> None:
        self.assertEqual(
            database_sources("kegg", "2026-02-01"),
            [
                "https://www.genome.jp/ftp/db/kofam/archives/2026-02-01/profiles.tar.gz",
                "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=",
            ],
        )

    def test_kegg_database_source_version_label_uses_requested_version(self) -> None:
        self.assertEqual(database_source_version_label("kegg", "2026-02-01"), "kofam archive 2026-02-01")

    def test_pfam_database_sources_use_requested_release_version(self) -> None:
        self.assertEqual(
            database_sources("pfam", "Pfam37.4"),
            [
                "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.4/Pfam-A.hmm.gz",
                "https://ecdm.loria.fr/data/EC-Pfam_calculated_associations_Extended.csv",
            ],
        )

    def test_pfam_database_source_version_label_uses_requested_version(self) -> None:
        self.assertEqual(database_source_version_label("pfam", "Pfam37.4"), "Pfam release Pfam37.4")

    def test_vfdb_database_source_version_label_uses_download_date(self) -> None:
        self.assertEqual(database_source_version_label("vfdb", "2026-04-24"), "VFDB_setB downloaded 2026-04-24")

    def test_amr_database_sources_use_requested_release_version(self) -> None:
        self.assertEqual(
            database_sources("amr", "2025-07-16.1"),
            [
                "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/2025-07-16.1/NCBIfam-AMRFinder.HMM.tar.gz",
                "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/2025-07-16.1/NCBIfam-AMRFinder.tsv",
            ],
        )

    def test_amr_database_source_version_label_uses_requested_version(self) -> None:
        self.assertEqual(database_source_version_label("amr", "2025-07-16.1"), "NCBIfam-AMRFinder 2025-07-16.1")

    def test_cazy_database_sources_use_requested_upstream_version(self) -> None:
        self.assertEqual(
            database_sources("cazy", "V14"),
            ["https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/V14/dbCAN-HMMdb-V14.txt"],
        )

    def test_cazy_database_source_version_label_uses_requested_version(self) -> None:
        self.assertEqual(database_source_version_label("cazy", "V14"), "dbCAN-HMMdb-V14")

    def test_replace_config_value_updates_single_key(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text('GTDB_DB: "/old/path"  # keep comment\nKEGG_DB: "/keep"\n', encoding="utf-8")

            from drakkar import cli as cli_module

            original = cli_module.CONFIG_PATH
            try:
                cli_module.CONFIG_PATH = config_path
                replace_config_value("GTDB_DB", "/new/path")
            finally:
                cli_module.CONFIG_PATH = original

            self.assertEqual(
                config_path.read_text(encoding="utf-8"),
                'GTDB_DB: "/new/path"  # keep comment\nKEGG_DB: "/keep"\n',
            )

    def test_set_default_database_path_writes_release_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text('KEGG_DB: "/old/path"\n', encoding="utf-8")

            from drakkar import cli as cli_module

            original = cli_module.CONFIG_PATH
            try:
                cli_module.CONFIG_PATH = config_path
                cli_module.set_default_database_path("kegg", "/tmp/kofams", "2026-02-01")
            finally:
                cli_module.CONFIG_PATH = original

            self.assertEqual(
                config_path.read_text(encoding="utf-8"),
                'KEGG_DB: "/tmp/kofams/2026-02-01"\n',
            )


if __name__ == "__main__":
    unittest.main()
