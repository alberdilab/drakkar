from __future__ import annotations

import subprocess
import sys
import unittest
from unittest.mock import patch

from drakkar import cli as cli_module


class UpdateCommandTests(unittest.TestCase):
    def test_run_update_displays_success_banner_with_installed_version(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=0)
        with patch.object(cli_module.subprocess, "run", return_value=completed) as run_mock:
            with patch.object(cli_module, "get_installed_drakkar_version", return_value="1.4.0"):
                with patch.object(cli_module, "display_update_success") as display_mock:
                    exit_code = cli_module.run_update()

        self.assertEqual(exit_code, 0)
        run_mock.assert_called_once()
        self.assertEqual(
            run_mock.call_args.args[0],
            [
                sys.executable, "-m", "pip", "install",
                "--upgrade", "--force-reinstall",
                "git+https://github.com/alberdilab/drakkar.git",
            ],
        )
        display_mock.assert_called_once_with("1.4.0")

    def test_run_update_skip_deps_uses_no_deps_flag(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=0)
        with patch.object(cli_module.subprocess, "run", return_value=completed) as run_mock:
            with patch.object(cli_module, "get_installed_drakkar_version", return_value="1.5.4"):
                with patch.object(cli_module, "display_update_success") as display_mock:
                    exit_code = cli_module.run_update(skip_deps=True)

        self.assertEqual(exit_code, 0)
        self.assertEqual(
            run_mock.call_args.args[0],
            [
                sys.executable, "-m", "pip", "install",
                "--upgrade", "--force-reinstall",
                "--no-deps",
                "git+https://github.com/alberdilab/drakkar.git",
            ],
        )
        display_mock.assert_called_once_with("1.5.4")

    def test_run_update_returns_failure_code_without_success_banner(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=3)
        with patch.object(cli_module.subprocess, "run", return_value=completed):
            with patch.object(cli_module, "display_update_success") as display_mock:
                exit_code = cli_module.run_update()

        self.assertEqual(exit_code, 3)
        display_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
