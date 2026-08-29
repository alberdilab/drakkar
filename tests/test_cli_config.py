from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from drakkar import cli as cli_module
from drakkar import config_persistence as persistence


class ConfigCommandTests(unittest.TestCase):
    def test_view_config_prints_path_and_contents(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text("TEST_KEY: value\n", encoding="utf-8")

            original = cli_module.CONFIG_PATH
            try:
                cli_module.CONFIG_PATH = config_path
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer):
                    exit_code = cli_module.view_config()
            finally:
                cli_module.CONFIG_PATH = original

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn(str(config_path.resolve()), output)
            self.assertIn("TEST_KEY: value", output)

    def test_resolve_editor_command_prefers_visual(self) -> None:
        with patch.dict(cli_module.os.environ, {"VISUAL": "code --wait"}, clear=False):
            self.assertEqual(cli_module.resolve_editor_command(), ["code", "--wait"])

    def test_edit_config_invokes_editor_with_config_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text("TEST_KEY: value\n", encoding="utf-8")

            original = cli_module.CONFIG_PATH
            try:
                cli_module.CONFIG_PATH = config_path
                with patch.object(cli_module, "resolve_editor_command", return_value=["/usr/bin/vi"]):
                    with patch.object(cli_module.subprocess, "run") as run_mock:
                        exit_code = cli_module.edit_config()
            finally:
                cli_module.CONFIG_PATH = original

            self.assertEqual(exit_code, 0)
            run_mock.assert_called_once_with(["/usr/bin/vi", str(config_path)], check=True)


class RestoreConfigTests(unittest.TestCase):
    @contextlib.contextmanager
    def _installation(self, config_text):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(config_text, encoding="utf-8")
            original = cli_module.CONFIG_PATH
            try:
                cli_module.CONFIG_PATH = config_path
                with patch.dict("os.environ", {"DRAKKAR_HOME": str(Path(tmpdir) / "home")}, clear=False):
                    yield config_path
            finally:
                cli_module.CONFIG_PATH = original

    def test_restore_writes_the_saved_paths_back_into_the_shipped_config(self) -> None:
        shipped = 'PFAM_DB: "/shipped/pfam/Pfam37.4"\nDREP_ANI: 0.98\n'
        with self._installation(shipped) as config_path:
            persistence.record_values({"PFAM_DB": "/site/pfam/Pfam38.2"})
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.restore_config()

            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 0)
        self.assertEqual(config["PFAM_DB"], "/site/pfam/Pfam38.2")
        self.assertEqual(config["DREP_ANI"], 0.98)
        self.assertIn("PFAM_DB", buffer.getvalue())

    def test_restore_reports_when_nothing_was_ever_saved(self) -> None:
        with self._installation('PFAM_DB: "/shipped/pfam/Pfam37.4"\n') as config_path:
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.restore_config()
            self.assertEqual(config_path.read_text(encoding="utf-8"), 'PFAM_DB: "/shipped/pfam/Pfam37.4"\n')

        self.assertEqual(exit_code, 1)
        self.assertIn("No saved configuration values", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
