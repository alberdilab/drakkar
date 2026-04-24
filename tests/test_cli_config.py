from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module


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


if __name__ == "__main__":
    unittest.main()
