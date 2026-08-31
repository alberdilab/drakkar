from __future__ import annotations

import re
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
        self.assertIn(
            '--default-resources "mem_mb=min(1048576, 24576 * 2 ** (attempt - 1))" '
            '"runtime=min(20160, 40 * 2 ** (attempt - 1))"',
            command,
        )

    def test_snakemake_runners_pass_resource_multipliers(self) -> None:
        cases = [
            lambda: cli_module.run_snakemake_environments("environments", "/tmp/envs", "local", 3, 4),
            lambda: cli_module.run_snakemake_preprocessing(
                "preprocessing", "project", "/tmp/output", False, "/tmp/envs", "local", False, False, False, 3, 4
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
                '--default-resources "mem_mb=min(1024, 24576 * 2 ** (attempt - 1))" '
                '"runtime=min(14, 40 * 2 ** (attempt - 1))" ',
            )

    def test_default_resources_escalate_across_attempts(self) -> None:
        """The defaults must grow on retry, and stay within the configured caps."""
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MAX_GB": 64, "SNAKEMAKE_MAX_TIME": 600}):
            args = cli_module.default_resource_args()
        exprs = dict(re.findall(r'"(\w+)=(.+?)"', args))
        for name, first, second in (("mem_mb", 8 * 1024, 16 * 1024), ("runtime", 10, 20)):
            self.assertEqual(eval(exprs[name], {"min": min, "attempt": 1}), first)
            self.assertEqual(eval(exprs[name], {"min": min, "attempt": 2}), second)
        # The cap still binds once doubling would exceed it.
        self.assertEqual(eval(exprs["runtime"], {"min": min, "attempt": 9}), 600)
        self.assertEqual(eval(exprs["mem_mb"], {"min": min, "attempt": 9}), 64 * 1024)

    def test_run_snakemake_cataloging_passes_binner_config(self) -> None:
        command = self.command_from_run(
            lambda: cli_module.run_snakemake_cataloging(
                "cataloging",
                "project",
                "/tmp/output",
                "/tmp/envs",
                "local",
                binners="metabat,comebin",
            )
        )

        self.assertIn("binners=metabat,comebin", command)

    def test_normalize_cataloging_binners(self) -> None:
        self.assertEqual(
            cli_module.normalize_cataloging_binners("semibin2,metabat,maxbin2"),
            "metabat,maxbin,semibin",
        )
        self.assertEqual(
            cli_module.normalize_cataloging_binners("all"),
            cli_module.DEFAULT_CATALOGING_BINNERS,
        )
        with patch.object(cli_module, "print"):
            self.assertIsNone(cli_module.normalize_cataloging_binners("metabat,unknown"))


if __name__ == "__main__":
    unittest.main()
