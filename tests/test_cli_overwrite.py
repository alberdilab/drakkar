from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module


class OverwriteCommandTests(unittest.TestCase):
    def test_prepare_output_directory_returns_true_when_not_locked(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.object(cli_module, "is_snakemake_locked", return_value=False):
                self.assertTrue(cli_module.prepare_output_directory(tmpdir, overwrite=False))

    def test_prepare_output_directory_overwrites_locked_directory_with_flag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "locked_output"
            output_path.mkdir()
            (output_path / "data.txt").write_text("content", encoding="utf-8")

            with patch.object(cli_module, "is_snakemake_locked", return_value=True):
                self.assertTrue(cli_module.prepare_output_directory(output_path, overwrite=True))

            self.assertFalse(output_path.exists())

    def test_prepare_output_directory_prompts_and_overwrites_locked_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "locked_output"
            output_path.mkdir()
            (output_path / "data.txt").write_text("content", encoding="utf-8")

            with patch.object(cli_module, "is_snakemake_locked", return_value=True):
                with patch.object(cli_module.sys.stdin, "isatty", return_value=True):
                    with patch("builtins.input", return_value="y"):
                        self.assertTrue(cli_module.prepare_output_directory(output_path, overwrite=False))

            self.assertFalse(output_path.exists())

    def test_prepare_output_directory_fails_non_interactively_without_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "locked_output"
            output_path.mkdir()

            buffer = io.StringIO()
            with patch.object(cli_module, "is_snakemake_locked", return_value=True):
                with patch.object(cli_module.sys.stdin, "isatty", return_value=False):
                    with contextlib.redirect_stdout(buffer):
                        result = cli_module.prepare_output_directory(output_path, overwrite=False)

            self.assertFalse(result)
            self.assertTrue(output_path.exists())
            self.assertIn("--overwrite", buffer.getvalue())
            self.assertIn("drakkar logging", buffer.getvalue())

    def test_validate_launch_metadata_directory_rejects_unwritable_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir)

            buffer = io.StringIO()
            with patch.object(cli_module.os, "access", return_value=False):
                with contextlib.redirect_stdout(buffer):
                    result = cli_module.validate_launch_metadata_directory(output_path)

            self.assertFalse(result)
            self.assertIn("Cannot write Drakkar run metadata", buffer.getvalue())
            self.assertIn("-o/--output", buffer.getvalue())

    def test_write_launch_metadata_returns_false_on_permission_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = type("Args", (), {"command": "annotating"})()

            buffer = io.StringIO()
            with patch.object(cli_module.os, "access", return_value=True):
                with patch("builtins.open", side_effect=PermissionError("denied")):
                    with contextlib.redirect_stdout(buffer):
                        result = cli_module.write_launch_metadata(args, tmpdir)

            self.assertFalse(result)
            self.assertIn("Cannot write Drakkar run metadata", buffer.getvalue())
            self.assertIn("PermissionError", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
