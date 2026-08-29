from __future__ import annotations

import contextlib
import subprocess
import sys
import unittest
from unittest.mock import patch

from drakkar import cli as cli_module


@contextlib.contextmanager
def isolated_config(saved=None, resolutions=()):
    """Keep the real config.yaml and ~/.drakkar out of the update tests."""
    with patch.object(cli_module, "snapshot_config", return_value=dict(saved or {})) as snapshot_mock:
        with patch.object(cli_module, "reconcile_config", return_value=list(resolutions)) as reconcile_mock:
            with patch.object(cli_module, "print_reconcile_report") as report_mock:
                yield snapshot_mock, reconcile_mock, report_mock


class UpdateCommandTests(unittest.TestCase):
    def test_run_update_displays_success_banner_with_installed_version(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=0)
        with isolated_config():
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
        with isolated_config():
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
        with isolated_config() as (_, reconcile_mock, _report):
            with patch.object(cli_module.subprocess, "run", return_value=completed):
                with patch.object(cli_module, "display_update_success") as display_mock:
                    exit_code = cli_module.run_update()

        self.assertEqual(exit_code, 3)
        display_mock.assert_not_called()
        # A failed install leaves the old config in place, so nothing is rewritten.
        reconcile_mock.assert_not_called()

    def test_run_update_saves_config_values_before_pip_and_restores_them_after(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=0)
        saved = {"PFAM_DB": "/site/pfam/Pfam38.2"}
        calls = []
        with isolated_config(saved=saved) as (snapshot_mock, reconcile_mock, report_mock):
            snapshot_mock.side_effect = lambda *args, **kwargs: calls.append("snapshot") or dict(saved)
            reconcile_mock.side_effect = lambda *args, **kwargs: calls.append("reconcile") or []
            with patch.object(cli_module.subprocess, "run", side_effect=lambda *a, **k: calls.append("pip") or completed):
                with patch.object(cli_module, "get_installed_drakkar_version", return_value="2.4.2"):
                    with patch.object(cli_module, "display_update_success"):
                        exit_code = cli_module.run_update()

        self.assertEqual(exit_code, 0)
        self.assertEqual(calls, ["snapshot", "pip", "reconcile"])
        self.assertEqual(reconcile_mock.call_args.kwargs, {"saved": saved})
        report_mock.assert_called_once()

    def test_run_update_still_installs_when_the_config_cannot_be_saved(self) -> None:
        completed = subprocess.CompletedProcess(args=["pip"], returncode=0)
        with isolated_config() as (snapshot_mock, reconcile_mock, _report):
            snapshot_mock.side_effect = OSError("read-only home")
            with patch.object(cli_module.subprocess, "run", return_value=completed):
                with patch.object(cli_module, "get_installed_drakkar_version", return_value="2.4.2"):
                    with patch.object(cli_module, "display_update_success") as display_mock:
                        exit_code = cli_module.run_update()

        self.assertEqual(exit_code, 0)
        display_mock.assert_called_once_with("2.4.2")
        reconcile_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
