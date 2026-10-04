from __future__ import annotations

import unittest
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.cli_parser import build_parser


class SkipDerepTests(unittest.TestCase):
    def command(self, **kwargs) -> str:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_profiling(
                    "profiling", "project", "genomes", "/tmp/output", "/tmp/envs", "local",
                    False, 0.98, False, None, **kwargs,
                )
        return run_mock.call_args.args[0][2]

    def test_profiling_dereplicates_by_default(self) -> None:
        self.assertIn("SKIP_DEREP=False", self.command())

    def test_skip_derep_reaches_the_workflow(self) -> None:
        self.assertIn("SKIP_DEREP=True", self.command(skip_derep=True))

    def test_the_profiling_parser_takes_skip_derep(self) -> None:
        args = build_parser().parse_args(["profiling", "-B", "bins.txt", "-R", "reads.tsv", "--skip-derep"])
        self.assertTrue(args.skip_derep)
        args = build_parser().parse_args(["profiling", "-B", "bins.txt", "-R", "reads.tsv"])
        self.assertFalse(args.skip_derep)


if __name__ == "__main__":
    unittest.main()
