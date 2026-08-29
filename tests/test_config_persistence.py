from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from drakkar import config_persistence as persistence

SHIPPED_CONFIG = """\
# Software
SNAKEMAKE_MODULE: "snakemake/8.25.3"

# Annotation databases
DATABASES_DIR: "/shipped/databases"
PFAM_DB: "/shipped/databases/pfam/Pfam37.4"
KEGG_DB: "/shipped/databases/kofams/2026-02-01"
AMRFINDER_DB: ""
GTDB_DB_226: "/shipped/gtdb/release226"

# Parameters
DREP_ANI: 0.98
"""


def make_release(directory, name):
    release = Path(directory) / name
    release.mkdir(parents=True)
    (release / "database_versions.yaml").write_text("requested_version: {}\n".format(name), encoding="utf-8")
    return release


class ConfigPersistenceTestCase(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self._tmpdir.name)
        self.home = self.root / "home"
        self.config_path = self.root / "config.yaml"
        self.config_path.write_text(SHIPPED_CONFIG, encoding="utf-8")
        patcher = patch.dict("os.environ", {"DRAKKAR_HOME": str(self.home)}, clear=False)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(self._tmpdir.cleanup)


class PreservedKeyTests(ConfigPersistenceTestCase):
    def test_database_and_directory_keys_are_preserved(self) -> None:
        for key in ("DATABASES_DIR", "ENVIRONMENTS_DIR", "PFAM_DB", "GTDB_DB_226", "PROSTT5_MODEL"):
            self.assertTrue(persistence.is_preserved_key(key), key)

    def test_software_and_parameter_keys_are_not_preserved(self) -> None:
        for key in ("SNAKEMAKE_MODULE", "DREP_ANI", "GENOMAD_SPLITS", "MEMORY_MULTIPLIER"):
            self.assertFalse(persistence.is_preserved_key(key), key)


class StoreTests(ConfigPersistenceTestCase):
    def test_record_values_round_trips_through_the_store(self) -> None:
        persistence.record_values({"PFAM_DB": "/site/pfam/Pfam38.2", "DREP_ANI": 0.98, "KEGG_DB": ""})
        self.assertEqual(persistence.read_store(), {"PFAM_DB": "/site/pfam/Pfam38.2"})
        self.assertTrue(persistence.store_path().is_file())

    def test_record_values_merges_into_earlier_values(self) -> None:
        persistence.record_values({"PFAM_DB": "/site/pfam/Pfam38.2"})
        persistence.record_values({"KEGG_DB": "/site/kofams/2026-07-02"})
        self.assertEqual(
            persistence.read_store(),
            {"PFAM_DB": "/site/pfam/Pfam38.2", "KEGG_DB": "/site/kofams/2026-07-02"},
        )

    def test_snapshot_prefers_the_current_config_over_earlier_saved_values(self) -> None:
        persistence.record_values({"PFAM_DB": "/stale/pfam/Pfam30.0", "VFDB_DB": "/site/vfdb/2026-04-24"})
        saved = persistence.snapshot_config(self.config_path)
        self.assertEqual(saved["PFAM_DB"], "/shipped/databases/pfam/Pfam37.4")
        # Keys the config no longer mentions stay in the store.
        self.assertEqual(saved["VFDB_DB"], "/site/vfdb/2026-04-24")

    def test_snapshot_writes_a_backup_of_the_whole_config(self) -> None:
        persistence.snapshot_config(self.config_path)
        backups = list(persistence.backup_dir().glob("config-*.yaml"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_text(encoding="utf-8"), SHIPPED_CONFIG)


class ResolveValueTests(ConfigPersistenceTestCase):
    def test_saved_path_wins_when_it_exists(self) -> None:
        databases = self.root / "site" / "pfam"
        saved = make_release(databases, "Pfam37.4")
        value, source, _ = persistence.resolve_value("PFAM_DB", str(saved), "/shipped/pfam/Pfam37.4")
        self.assertEqual(value, str(saved))
        self.assertEqual(source, persistence.SOURCE_CONFIG)

    def test_newer_release_beside_the_saved_one_wins(self) -> None:
        databases = self.root / "site" / "pfam"
        saved = make_release(databases, "Pfam37.4")
        newest = make_release(databases, "Pfam38.2")
        value, source, _ = persistence.resolve_value("PFAM_DB", str(saved), "/shipped/pfam/Pfam37.4")
        self.assertEqual(value, str(newest))
        self.assertEqual(source, persistence.SOURCE_NEWEST)

    def test_missing_saved_path_falls_back_to_the_newest_release_that_exists(self) -> None:
        databases = self.root / "site" / "kofams"
        make_release(databases, "2026-02-01")
        newest = make_release(databases, "2026-07-02")
        value, source, _ = persistence.resolve_value(
            "KEGG_DB", str(databases / "2025-01-01"), "/shipped/kofams/2026-02-01"
        )
        self.assertEqual(value, str(newest))
        self.assertEqual(source, persistence.SOURCE_NEWEST)

    def test_directories_without_a_release_marker_are_ignored(self) -> None:
        databases = self.root / "site" / "pfam"
        saved = make_release(databases, "Pfam37.4")
        (databases / "tmp-9999").mkdir()
        value, _, _ = persistence.resolve_value("PFAM_DB", str(saved), "/shipped/pfam/Pfam37.4")
        self.assertEqual(value, str(saved))

    def test_saved_value_is_kept_when_nothing_exists_on_disk(self) -> None:
        value, source, detail = persistence.resolve_value(
            "PFAM_DB", "/unmounted/pfam/Pfam38.2", "/shipped/pfam/Pfam37.4"
        )
        self.assertEqual(value, "/unmounted/pfam/Pfam38.2")
        self.assertEqual(source, persistence.SOURCE_CONFIG)
        self.assertIn("no release found on disk", detail)

    def test_pinned_release_keys_are_never_repointed(self) -> None:
        databases = self.root / "site" / "gtdb"
        make_release(databases, "release226")
        make_release(databases, "release232")
        value, source, _ = persistence.resolve_value(
            "GTDB_DB_226", str(databases / "release226"), "/shipped/gtdb/release226"
        )
        self.assertEqual(value, str(databases / "release226"))
        self.assertEqual(source, persistence.SOURCE_CONFIG)

    def test_shipped_default_is_used_when_nothing_was_saved(self) -> None:
        value, source, _ = persistence.resolve_value("PFAM_DB", "", "/shipped/pfam/Pfam37.4")
        self.assertEqual(value, "/shipped/pfam/Pfam37.4")
        self.assertEqual(source, persistence.SOURCE_INSTALLED)


class ReconcileTests(ConfigPersistenceTestCase):
    def test_saved_values_are_written_back_into_a_replaced_config(self) -> None:
        databases = self.root / "site"
        pfam = make_release(databases / "pfam", "Pfam38.2")
        persistence.record_values(
            {
                "DATABASES_DIR": str(databases),
                "PFAM_DB": str(pfam),
                "KEGG_DB": "/unmounted/kofams/2026-07-02",
                "GTDB_DB_226": "/site/gtdb/release226",
            }
        )

        resolutions = persistence.reconcile_config(self.config_path)

        config = yaml.safe_load(self.config_path.read_text(encoding="utf-8"))
        self.assertEqual(config["DATABASES_DIR"], str(databases))
        self.assertEqual(config["PFAM_DB"], str(pfam))
        self.assertEqual(config["KEGG_DB"], "/unmounted/kofams/2026-07-02")
        self.assertEqual(config["GTDB_DB_226"], "/site/gtdb/release226")
        # Untouched keys keep the values, comments and layout of the new release.
        self.assertEqual(config["SNAKEMAKE_MODULE"], "snakemake/8.25.3")
        self.assertEqual(config["DREP_ANI"], 0.98)
        self.assertIn("# Annotation databases", self.config_path.read_text(encoding="utf-8"))
        self.assertEqual({item.key for item in resolutions if item.changed}, {"DATABASES_DIR", "PFAM_DB", "KEGG_DB", "GTDB_DB_226"})

    def test_keys_dropped_by_the_new_release_are_not_reintroduced(self) -> None:
        persistence.record_values({"FOLDSEEK_DB": "/site/foldseek/2026-04-24/foldseek_db"})
        persistence.reconcile_config(self.config_path)
        self.assertNotIn("FOLDSEEK_DB", self.config_path.read_text(encoding="utf-8"))

    def test_reconciliation_saves_the_values_it_settled_on(self) -> None:
        persistence.reconcile_config(self.config_path)
        self.assertEqual(persistence.read_store()["PFAM_DB"], "/shipped/databases/pfam/Pfam37.4")

    def test_report_lists_restored_values(self) -> None:
        persistence.record_values({"PFAM_DB": "/site/pfam/Pfam38.2"})
        resolutions = persistence.reconcile_config(self.config_path)
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            persistence.print_reconcile_report(resolutions, self.config_path)
        output = buffer.getvalue()
        self.assertIn("PFAM_DB: /site/pfam/Pfam38.2", output)
        self.assertIn(str(persistence.store_path()), output)


if __name__ == "__main__":
    unittest.main()
