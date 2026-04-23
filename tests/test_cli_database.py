from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from drakkar.cli import replace_config_value, validate_database_version
from drakkar.database_registry import (
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

    def test_database_release_dir_joins_base_and_version(self) -> None:
        release_dir = database_release_dir("amr", "/tmp/amr", "20260421")
        self.assertEqual(release_dir, Path("/tmp/amr/20260421"))

    def test_database_target_path_uses_database_specific_basename(self) -> None:
        target_path = database_target_path("kegg", "/tmp/kegg", "20260421")
        self.assertEqual(target_path, Path("/tmp/kegg/20260421/kofams"))

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


if __name__ == "__main__":
    unittest.main()
