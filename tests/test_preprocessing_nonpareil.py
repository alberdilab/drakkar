from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from drakkar import cli as cli_module


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
NONPAREIL_RULES = ROOT / "drakkar" / "workflow" / "rules" / "preprocessing_nonpareil.smk"
PREPROCESSING_SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "preprocessing_stats.py"


class PreprocessingNonpareilTests(unittest.TestCase):
    def test_run_snakemake_preprocessing_passes_nonpareil_config(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module.subprocess, "run") as run_mock:
                cli_module.run_snakemake_preprocessing(
                    "preprocessing",
                    "project",
                    "/tmp/output",
                    False,
                    "/tmp/envs",
                    "local",
                    nonpareil=True,
                )

        command = run_mock.call_args.args[0][2]
        self.assertIn("nonpareil=True", command)

    def test_run_snakemake_preprocessing_defaults_nonpareil_to_false(self) -> None:
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
        self.assertIn("nonpareil=False", command)

    def test_preprocessing_targets_nonpareil_outputs(self) -> None:
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        self.assertIn('NONPAREIL = config.get("nonpareil", None)', snakefile)
        self.assertGreaterEqual(snakefile.count('if NONPAREIL else []'), 2)
        self.assertIn("rules/preprocessing_nonpareil.smk", snakefile)
        self.assertIn("preprocessing/nonpareil/{{sample}}_np.tsv", snakefile)

    def test_preprocessing_has_nonpareil_rule(self) -> None:
        rules = NONPAREIL_RULES.read_text(encoding="utf-8")
        self.assertIn("rule nonpareil:", rules)
        self.assertIn("NONPAREIL_MODULE", rules)
        self.assertIn("nonpareil -s", rules)
        self.assertIn("nonpareil_stats.R", rules)
        self.assertIn("_np.tsv", rules)

    def test_preprocessing_stats_table_includes_requested_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fastp_json = tmp_path / "sample1.json"
            fastp_json.write_text(
                json.dumps(
                    {
                        "summary": {
                            "before_filtering": {"total_reads": 100, "total_bases": 1000},
                            "after_filtering": {"total_reads": 80, "total_bases": 800},
                        },
                        "adapter_cutting": {
                            "adapter_trimmed_reads": 10,
                            "adapter_trimmed_bases": 50,
                        },
                    }
                ),
                encoding="utf-8",
            )

            fastq = tmp_path / "sample1_1.fq.gz"
            with gzip.open(fastq, "wt", encoding="utf-8") as handle:
                handle.write("@r1\nACGT\n+\n!!!!\n@r2\nAC\n+\n!!\n")

            smf = tmp_path / "sample1_smf.tsv"
            smf.write_text("sample\tread_fraction\nsample1\t0.85\n", encoding="utf-8")

            nonpareil = tmp_path / "sample1_np.tsv"
            nonpareil.write_text(
                "sample\tkappa\tC\tLR\tmodelR\tLRstar\tdiversity\n"
                "sample1\t0.1\t0.91\t1.2\t0.98\t2.3\t4.5\n",
                encoding="utf-8",
            )

            output = tmp_path / "preprocessing.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(PREPROCESSING_SCRIPT),
                    "-p",
                    str(fastp_json),
                    "-f",
                    str(fastq),
                    "-s",
                    str(smf),
                    "-n",
                    str(nonpareil),
                    "-o",
                    str(output),
                ],
                check=True,
            )

            table = pd.read_csv(output, sep="\t")
            self.assertEqual(
                list(table.columns),
                [
                    "sample",
                    "reads_pre_fastp",
                    "bases_pre_fastp",
                    "adapter_trimmed_reads",
                    "adapter_trimmed_bases",
                    "reads_post_fastp",
                    "bases_post_fastp",
                    "host_reads",
                    "host_bases",
                    "metagenomic_reads",
                    "metagenomic_bases",
                    "singlem_fraction",
                    "nonpareil_C",
                    "nonpareil_LR",
                    "nonpareil_modelR",
                    "nonpareil_LRstar",
                    "nonpareil_diversity",
                ],
            )
            row = table.iloc[0].to_dict()
            self.assertEqual(row["sample"], "sample1")
            self.assertEqual(row["reads_pre_fastp"], 100)
            self.assertEqual(row["adapter_trimmed_reads"], 10)
            self.assertEqual(row["metagenomic_reads"], 4)
            self.assertEqual(row["metagenomic_bases"], 12)
            self.assertEqual(row["singlem_fraction"], 0.85)
            self.assertEqual(row["nonpareil_C"], 0.91)
            self.assertEqual(row["nonpareil_modelR"], 0.98)


if __name__ == "__main__":
    unittest.main()
