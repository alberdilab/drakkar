from __future__ import annotations

import contextlib
import io
import unittest
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.cli import available_gtdb_versions, validate_gtdb_version


class GtdbVersionTests(unittest.TestCase):
    def test_available_gtdb_versions_uses_numbered_config_keys(self) -> None:
        config = {
            "GTDB_DB": "/default",
            "GTDB_DB_220": "/release220",
            "GTDB_DB_232": "/release232",
            "GTDB_DB_DEV": "/ignored",
        }

        self.assertEqual(available_gtdb_versions(config), ["232", "220"])

    def test_validate_gtdb_version_accepts_configured_release(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(validate_gtdb_version("226", {"GTDB_DB_226": "/release226"}), "226")

    def test_validate_gtdb_version_rejects_unconfigured_release(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNone(validate_gtdb_version("999", {"GTDB_DB_226": "/release226"}))

    def test_run_snakemake_annotating_passes_requested_gtdb_version(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module.subprocess, "run") as run_mock:
                cli_module.run_snakemake_annotating(
                    "annotating",
                    "project",
                    "taxonomy",
                    "/tmp/output",
                    "/tmp/envs",
                    "local",
                    gtdb_version="226",
                )

        command = run_mock.call_args.args[0][2]
        self.assertIn("gtdb_version=226", command)

    def test_run_snakemake_annotating_uses_default_gtdb_when_version_omitted(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module.subprocess, "run") as run_mock:
                cli_module.run_snakemake_annotating(
                    "annotating",
                    "project",
                    "taxonomy",
                    "/tmp/output",
                    "/tmp/envs",
                    "local",
                )

        command = run_mock.call_args.args[0][2]
        self.assertNotIn("gtdb_version=", command)


if __name__ == "__main__":
    unittest.main()
