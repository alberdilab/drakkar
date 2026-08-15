from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from drakkar import cli as cli_module
from drakkar.cli_parser import build_parser
from drakkar.cli_validation import bin_filter_args


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "drakkar" / "workflow" / "config.yaml"
CATALOGING_RULES = ROOT / "drakkar" / "workflow" / "rules" / "cataloging.smk"
CATALOGING_ENV = ROOT / "drakkar" / "workflow" / "envs" / "cataloging.yaml"
SCRIPTS_DIR = ROOT / "drakkar" / "workflow" / "scripts"
ALL_BIN_PATHS_SCRIPT = SCRIPTS_DIR / "all_bin_paths.py"
ALL_BIN_METADATA_SCRIPT = SCRIPTS_DIR / "all_bin_metadata.py"

sys.path.insert(0, str(SCRIPTS_DIR))

from bin_report import bin_id_column, bin_ids_from_table  # noqa: E402

BINETTE_REPORT = (
    "name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tscore\t"
    "checkm2_model\tsize\tN50\tcoding_density\tcontig_count\n"
    "binette_bin1\tbinette\tFalse\tmetabat_3\t95.1\t2.0\t91.1\tgeneral\t3000000\t50000\t0.88\t120\n"
    "binette_bin2\tmaxbin\tTrue\tmaxbin_7\t78.4\t4.5\t69.4\tgeneral\t2100000\t22000\t0.9\t310\n"
)


class BinFilterConfigTests(unittest.TestCase):
    def test_binette_rule_passes_all_filtering_thresholds(self) -> None:
        rules = CATALOGING_RULES.read_text(encoding="utf-8")
        binette_rule = rules.split("checkpoint binette:", 1)[1].split("# Regenerate the bin_id wildcard", 1)[0]

        self.assertIn("--prefix {params.bin_prefix}", binette_rule)
        self.assertIn("--min_completeness {params.min_completeness}", binette_rule)
        self.assertIn("--max_contamination {params.max_contamination}", binette_rule)
        self.assertIn("--min_length {params.min_bin_length}", binette_rule)
        self.assertIn("--max_length {params.max_bin_length}", binette_rule)

    def test_config_defines_bin_length_thresholds(self) -> None:
        config = CONFIG_PATH.read_text(encoding="utf-8")

        self.assertIn("MIN_COMPLETENESS: 70", config)
        self.assertIn("MAX_CONTAMINATION: 10", config)
        self.assertIn("MIN_BIN_LENGTH: 200000", config)
        self.assertIn("MAX_BIN_LENGTH: 10000000", config)

    def test_cataloging_environment_pins_binette_with_filtering_options(self) -> None:
        self.assertIn("binette==1.2.1", CATALOGING_ENV.read_text(encoding="utf-8"))

    def test_bin_fasta_paths_use_the_binette_bin_prefix(self) -> None:
        rules = CATALOGING_RULES.read_text(encoding="utf-8")

        self.assertIn("final_bins/{BINETTE_BIN_PREFIX}_bin{wildcards.bin_id}.fa", rules)


class BinFilterArgumentTests(unittest.TestCase):
    def test_cataloging_accepts_bin_filtering_arguments(self) -> None:
        args = build_parser().parse_args(
            [
                "cataloging",
                "-i",
                "reads",
                "--min-completeness",
                "80",
                "--max-contamination",
                "5",
                "--min-length",
                "150000",
                "--max-length",
                "12000000",
            ]
        )

        self.assertEqual(args.min_completeness, 80)
        self.assertEqual(args.max_contamination, 5)
        self.assertEqual(args.min_bin_length, 150000)
        self.assertEqual(args.max_bin_length, 12000000)

    def test_bin_filtering_arguments_default_to_config_values(self) -> None:
        args = build_parser().parse_args(["cataloging", "-i", "reads"])

        self.assertIsNone(args.min_completeness)
        self.assertIsNone(args.max_contamination)
        self.assertIsNone(args.min_bin_length)
        self.assertIsNone(args.max_bin_length)
        self.assertEqual(bin_filter_args(), "")

    def test_bin_filter_args_only_overrides_provided_values(self) -> None:
        overrides = bin_filter_args(min_completeness=80, max_bin_length=12000000)

        self.assertEqual(overrides, "MIN_COMPLETENESS=80 MAX_BIN_LENGTH=12000000 ")

    def test_completeness_and_contamination_must_be_percentages(self) -> None:
        parser = build_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args(["cataloging", "-i", "reads", "--max-contamination", "120"])


class BinFilterLauncherTests(unittest.TestCase):
    def test_cataloging_run_forwards_requested_thresholds(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_cataloging(
                    "cataloging",
                    "project",
                    "/tmp/output",
                    "/tmp/envs",
                    "local",
                    min_completeness=80,
                    max_contamination=5,
                    min_bin_length=150000,
                    max_bin_length=12000000,
                )

        command = run_mock.call_args.args[0][2]
        self.assertIn("MIN_COMPLETENESS=80", command)
        self.assertIn("MAX_CONTAMINATION=5", command)
        self.assertIn("MIN_BIN_LENGTH=150000", command)
        self.assertIn("MAX_BIN_LENGTH=12000000", command)

    def test_cataloging_run_omits_thresholds_left_at_config_defaults(self) -> None:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_cataloging(
                    "cataloging",
                    "project",
                    "/tmp/output",
                    "/tmp/envs",
                    "local",
                )

        command = run_mock.call_args.args[0][2]
        self.assertNotIn("MIN_COMPLETENESS", command)
        self.assertNotIn("MAX_BIN_LENGTH", command)


class BinReportParsingTests(unittest.TestCase):
    def test_bin_ids_are_parsed_from_binette_names(self) -> None:
        table = pd.DataFrame({"name": ["binette_bin1", "binette_bin12"]})

        self.assertEqual(bin_ids_from_table(table), ["1", "12"])
        self.assertEqual(bin_id_column(table), ["1", "12"])

    def test_legacy_bin_id_reports_are_still_supported(self) -> None:
        table = pd.DataFrame({"bin_id": [1, 2]})

        self.assertEqual(bin_ids_from_table(table), ["1", "2"])

    def test_reports_without_bin_columns_yield_no_ids(self) -> None:
        self.assertEqual(bin_ids_from_table(pd.DataFrame({"completeness": [90]})), [])

    def test_bin_summary_reports_the_drakkar_bin_id(self) -> None:
        sys.path.insert(0, str(SCRIPTS_DIR))
        import cataloging_stats

        with tempfile.TemporaryDirectory() as tmpdir:
            report = Path(tmpdir) / "assembly1.tsv"
            report.write_text(BINETTE_REPORT, encoding="utf-8")

            assembly, summary = cataloging_stats.read_bin_summary(str(report))

            self.assertEqual(assembly, "assembly1")
            self.assertEqual(summary["final_bins"], 2)
            self.assertEqual(summary["best_bin"], "1")
            self.assertEqual(summary["best_bin_score"], 91.1)

    def test_bin_exports_use_binette_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            report = tmp_path / "assembly1.tsv"
            report.write_text(BINETTE_REPORT, encoding="utf-8")
            paths = tmp_path / "all_bin_paths.txt"
            metadata = tmp_path / "all_bin_metadata.csv"

            subprocess.run(
                [sys.executable, str(ALL_BIN_PATHS_SCRIPT), str(report), "-o", str(paths)],
                check=True,
            )
            subprocess.run(
                [sys.executable, str(ALL_BIN_METADATA_SCRIPT), str(report), "-o", str(metadata)],
                check=True,
            )

            self.assertEqual(
                paths.read_text(encoding="utf-8").split(),
                [
                    "cataloging/final/assembly1/assembly1_bin_1.fa",
                    "cataloging/final/assembly1/assembly1_bin_2.fa",
                ],
            )

            table = pd.read_csv(metadata)
            self.assertEqual(
                table["genome"].tolist(),
                ["assembly1_bin_1.fa", "assembly1_bin_2.fa"],
            )
            self.assertEqual(table["completeness"].tolist(), [95.1, 78.4])
            self.assertEqual(table["contig_count"].tolist(), [120, 310])


if __name__ == "__main__":
    unittest.main()
