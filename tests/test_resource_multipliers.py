from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module


class ResourceMultiplierTests(unittest.TestCase):
    def command_from_run(self, run_callable):
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                run_callable()
        return run_mock.call_args.args[0][2]

    def assert_resource_overrides(self, command: str) -> None:
        self.assertIn("memory_multiplier=3", command)
        self.assertIn("time_multiplier=4", command)
        self.assertIn("--default-resources mem_mb=24576 runtime=40", command)

    def test_snakemake_runners_pass_resource_multipliers(self) -> None:
        cases = [
            lambda: cli_module.run_snakemake_environments("environments", "/tmp/envs", "local", 3, 4),
            lambda: cli_module.run_snakemake_preprocessing(
                "preprocessing", "project", "/tmp/output", False, "/tmp/envs", "local", False, False, 3, 4
            ),
            lambda: cli_module.run_snakemake_cataloging("cataloging", "project", "/tmp/output", "/tmp/envs", "local", 3, 4),
            lambda: cli_module.run_snakemake_profiling(
                "profiling", "project", "genomes", "/tmp/output", "/tmp/envs", "local", False, 0.98, False, None, 3, 4
            ),
            lambda: cli_module.run_snakemake_dereplicating(
                "dereplicating", "project", "/tmp/output", "/tmp/envs", "local", 0.98, False, None, 3, 4
            ),
            lambda: cli_module.run_snakemake_annotating(
                "annotating", "project", "taxonomy", "/tmp/output", "/tmp/envs", "local", None, 3, 4
            ),
            lambda: cli_module.run_snakemake_inspecting("inspecting", "project", "/tmp/output", "/tmp/envs", "local", 3, 4),
            lambda: cli_module.run_snakemake_expressing("expressing", "project", "/tmp/output", "/tmp/envs", "local", 3, 4),
            lambda: cli_module.run_snakemake_database(
                "database", "project", "/tmp/output", "/tmp/envs", "local", "kegg", Path("/tmp/db"), "2026-02-01", 120, 3, 4
            ),
        ]

        for case in cases:
            with self.subTest(case=case):
                self.assert_resource_overrides(self.command_from_run(case))

    def test_default_resources_are_capped_by_config(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MAX_GB": 1, "SNAKEMAKE_MAX_TIME": 14}):
            self.assertEqual(
                cli_module.default_resource_args(memory_multiplier=3, time_multiplier=4),
                "--default-resources mem_mb=1024 runtime=14 ",
            )


if __name__ == "__main__":
    unittest.main()
