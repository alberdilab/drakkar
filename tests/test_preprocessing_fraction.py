from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
PREPROCESSING_RULES = ROOT / "drakkar" / "workflow" / "rules" / "preprocessing.smk"
PREPROCESSING_REF_RULES = ROOT / "drakkar" / "workflow" / "rules" / "preprocessing_ref.smk"


class PreprocessingFractionTests(unittest.TestCase):
    def test_run_snakemake_preprocessing_passes_fraction_config(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module.subprocess, "run") as run_mock:
                cli_module.run_snakemake_preprocessing(
                    "preprocessing",
                    "project",
                    "/tmp/output",
                    False,
                    "/tmp/envs",
                    "local",
                    fraction=True,
                )

        command = run_mock.call_args.args[0][2]
        self.assertIn("fraction=True", command)

    def test_run_snakemake_preprocessing_defaults_fraction_to_false(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module.subprocess, "run") as run_mock:
                cli_module.run_snakemake_preprocessing(
                    "preprocessing",
                    "project",
                    "/tmp/output",
                    False,
                    "/tmp/envs",
                    "local",
                )

        command = run_mock.call_args.args[0][2]
        self.assertIn("fraction=False", command)

    def test_preprocessing_targets_singlem_fraction_outputs(self) -> None:
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        self.assertGreaterEqual(snakefile.count('if FRACTION else []'), 2)
        self.assertIn('preprocessing/singlem/{{sample}}_smf.tsv', snakefile)

    def test_preprocessing_without_reference_has_singlem_rules(self) -> None:
        rules = PREPROCESSING_RULES.read_text(encoding="utf-8")
        self.assertIn("rule singlem:", rules)
        self.assertIn("rule singlem_mf:", rules)
        self.assertIn('SINGLEM_DB = config["SINGLEM_DB"]', rules)
        self.assertIn("singlem pipe", rules)
        self.assertIn("singlem microbial_fraction", rules)
        self.assertEqual(rules.count("--metapackage {params.singlem_db}"), 2)

    def test_preprocessing_with_reference_has_singlem_rules(self) -> None:
        rules = PREPROCESSING_REF_RULES.read_text(encoding="utf-8")
        self.assertIn("rule singlem:", rules)
        self.assertIn("rule singlem_mf:", rules)
        self.assertIn('SINGLEM_DB = config["SINGLEM_DB"]', rules)
        self.assertIn("singlem pipe", rules)
        self.assertIn("singlem microbial_fraction", rules)
        self.assertEqual(rules.count("--metapackage {params.singlem_db}"), 2)


if __name__ == "__main__":
    unittest.main()
