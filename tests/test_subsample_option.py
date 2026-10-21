"""CLI, fraction, and preprocessing wiring for paired read subsampling."""

import argparse
import csv
import json
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.cli_parser import build_parser
from drakkar.input_manifests import (
    argument_preprocessed_to_json,
    argument_transcriptome_to_json,
    file_preprocessed_to_json,
    file_samples_to_json,
    file_transcriptome_to_json,
)
from drakkar.subsample import parse_subsample_value, subsample_fraction
from drakkar.workflow.scripts import subsample_reads


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "drakkar" / "workflow"


class SubsampleValueTests(unittest.TestCase):
    def test_supported_forms(self):
        cases = {
            "0.25": ("fraction", 0.25),
            ".5": ("fraction", 0.5),
            "50000000": ("pairs", 50_000_000),
            "50K": ("pairs", 50_000),
            "50M": ("pairs", 50_000_000),
            "2G": ("pairs", 2_000_000_000),
            "10Mb": ("bases", 10_000_000),
            "10Gb": ("bases", 10_000_000_000),
            "1.5Gb": ("bases", 1_500_000_000),
        }
        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(parse_subsample_value(value), expected)

    def test_invalid_forms(self):
        for value in ("0", "0.0", "1.0", "2.0", "-1", "0M", "0Gb", "50T", "10GB", "10mb", "wat"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    parse_subsample_value(value)
        with self.assertRaisesRegex(argparse.ArgumentTypeError, "omit --subsample"):
            parse_subsample_value("1.0")

    def test_cli_accepts_both_commands(self):
        for command in ("preprocessing", "complete", "cataloging", "profiling", "expressing"):
            for value in ("0.25", "50M", "10Gb"):
                with self.subTest(command=command, value=value):
                    table_arg = "-R" if command in ("profiling", "expressing") else "-f"
                    args = build_parser().parse_args([command, table_arg, "samples.tsv", "--subsample", value])
                    self.assertEqual(args.subsample, value)
            with self.assertRaises(SystemExit):
                table_arg = "-R" if command in ("profiling", "expressing") else "-f"
                build_parser().parse_args([command, table_arg, "samples.tsv", "--subsample", "1.0"])

    def test_fraction_never_upsamples(self):
        self.assertEqual(subsample_fraction("0.25", 10, 2000), 0.25)
        self.assertEqual(subsample_fraction("50M", 10, 2000), 1)
        self.assertEqual(subsample_fraction("10Gb", 10, 2000), 1)
        self.assertEqual(subsample_fraction("50K", 100_000, 2000), 0.5)
        self.assertEqual(subsample_fraction("10Mb", 100_000, 20_000_000), 0.5)
        self.assertEqual(subsample_fraction("50M", 0, 0), 1)

    def test_target_reaches_snakemake_config(self):
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_preprocessing(
                    "preprocessing", "project", "/tmp/output", False, "/tmp/envs", "local",
                    subsample="10Gb",
                )
        self.assertIn("subsample=10Gb", shlex.split(run_mock.call_args.args[0][2]))
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_profiling(
                    "profiling", "project", "genomes", "/tmp/output", "/tmp/envs", "local",
                    False, 0.98, False, None, subsample="50M",
                )
        self.assertIn("subsample=50M", shlex.split(run_mock.call_args.args[0][2]))
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_cataloging(
                    "cataloging", "project", "/tmp/output", "/tmp/envs", "local", subsample="10Gb",
                )
        self.assertIn("subsample=10Gb", shlex.split(run_mock.call_args.args[0][2]))
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_expressing(
                    "expressing", "project", "/tmp/output", "/tmp/envs", "local", subsample="0.25",
                )
        self.assertIn("subsample=0.25", shlex.split(run_mock.call_args.args[0][2]))

    def test_table_targets_override_fallback_and_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for sample in ("S1", "S2", "S3"):
                for mate in (1, 2):
                    (root / f"{sample}_{mate}.fq.gz").write_bytes(b"reads")
            table = root / "samples.tsv"
            table.write_text("sample\trawreads1\trawreads2\tsubsample\n" + "".join(
                f"{s}\t{root / (s + '_1.fq.gz')}\t{root / (s + '_2.fq.gz')}\t{v}\n"
                for s, v in (("S1", "10Gb"), ("S2", "NA"), ("S3", ""))
            ))
            file_samples_to_json(str(table), str(root), subsample="0.25")
            self.assertEqual(json.loads((root / "data/sample_to_subsample.json").read_text()),
                             {"S1": "10Gb", "S2": "0.25", "S3": "0.25"})
            file_samples_to_json(str(table), str(root))
            self.assertEqual(json.loads((root / "data/sample_to_subsample.json").read_text()), {"S1": "10Gb"})
            # Plain integers next to blank cells must stay text, not become 50000000.0.
            table.write_text(table.read_text().replace("\t10Gb\n", "\t50000000\n"))
            file_samples_to_json(str(table), str(root))
            self.assertEqual(json.loads((root / "data/sample_to_subsample.json").read_text()), {"S1": "50000000"})
            table.write_text(table.read_text().replace("\t50000000\n", "\t10Gb\n"))
            file_preprocessed_to_json(str(table), str(root), subsample="50M", enable_subsample=True)
            self.assertEqual(json.loads((root / "data/sample_to_subsample.json").read_text()),
                             {"S1": "10Gb", "S2": "50M", "S3": "50M"})
            file_preprocessed_to_json(str(table), str(root))
            self.assertEqual(json.loads((root / "data/sample_to_subsample.json").read_text()), {})
            table.write_text(table.read_text().replace("\t10Gb\n", "\tbad-target\n"))
            with patch("drakkar.input_manifests.report_input_resolution_errors", side_effect=ValueError) as report:
                with self.assertRaises(ValueError):
                    file_samples_to_json(str(table), str(root))
            self.assertIn("S1", report.call_args.args[0][0])
            self.assertIn("bad-target", report.call_args.args[0][0])
            with patch("drakkar.input_manifests.report_input_resolution_errors", side_effect=ValueError) as report:
                with self.assertRaises(ValueError):
                    file_preprocessed_to_json(str(table), str(root), enable_subsample=True)
            self.assertIn("bad-target", report.call_args.args[0][0])

    def test_cataloging_and_expressing_manifests(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for sample in ("S1", "S2"):
                for mate in (1, 2):
                    (root / f"{sample}_{mate}.fq.gz").write_bytes(b"reads")
            table = root / "samples.tsv"
            table.write_text("sample\trawreads1\trawreads2\tsubsample\n" + "".join(
                f"{sample}\t{root / (sample + '_1.fq.gz')}\t{root / (sample + '_2.fq.gz')}\t{value}\n"
                for sample, value in (("S1", "10Gb"), ("S2", "NA"))
            ))
            file_preprocessed_to_json(str(table), str(root), subsample="0.25", enable_subsample=True)
            target_path = root / "data/sample_to_subsample.json"
            self.assertEqual(json.loads(target_path.read_text()), {"S1": "10Gb", "S2": "0.25"})
            # Cataloging -i plus -f reads files from the directory but keeps
            # per-sample target metadata from the table.
            argument_preprocessed_to_json(str(root), str(root), subsample="0.25", infofile=str(table))
            self.assertEqual(json.loads(target_path.read_text()), {"S1": "10Gb", "S2": "0.25"})
            argument_preprocessed_to_json(str(root), str(root))
            self.assertEqual(json.loads(target_path.read_text()), {})
            file_transcriptome_to_json(str(table), str(root), subsample="50M")
            self.assertEqual(json.loads(target_path.read_text()), {"S1": "10Gb", "S2": "50M"})
            file_transcriptome_to_json(str(table), str(root))
            self.assertEqual(json.loads(target_path.read_text()), {"S1": "10Gb"})
            argument_transcriptome_to_json(str(root), str(root), subsample="50M")
            self.assertEqual(json.loads(target_path.read_text()), {"S1": "50M", "S2": "50M"})
            file_preprocessed_to_json(str(table), str(root))
            self.assertEqual(json.loads(target_path.read_text()), {})
            table.write_text(table.read_text().replace("10Gb", "invalid"))
            with patch("drakkar.input_manifests.report_input_resolution_errors", side_effect=ValueError) as report:
                with self.assertRaises(ValueError):
                    file_transcriptome_to_json(str(table), str(root))
            self.assertIn("S1", report.call_args.args[0][0])
            self.assertIn("invalid", report.call_args.args[0][0])
            with patch("drakkar.input_manifests.report_input_resolution_errors", side_effect=ValueError) as report:
                with self.assertRaises(ValueError):
                    argument_preprocessed_to_json(str(root), str(root), infofile=str(table))
            self.assertIn("invalid", report.call_args.args[0][0])

    def test_summary_reads_are_both_mates_and_na_when_unused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastp = root / "S1.json"
            fastp.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 8, "total_bases": 800}, "after_filtering": {"total_reads": 6, "total_bases": 600}}}))
            stats = root / "S1_subsample.tsv"
            stats.write_text("sample\treads_before\treads_subsampled\tpairs_before\tpairs_after\nS1\t20\t8\t10\t4\n")
            output = root / "preprocessing.tsv"
            command = [sys.executable, str(WORKFLOW / "scripts" / "preprocessing_stats.py"), "-p", str(fastp), "-o", str(output)]
            subprocess.run(command + ["--subsample", str(stats)], check=True)
            with output.open() as handle:
                self.assertEqual(next(csv.DictReader(handle, delimiter="\t"))["reads_subsampled"], "8")
            subprocess.run(command, check=True)
            with output.open() as handle:
                self.assertEqual(next(csv.DictReader(handle, delimiter="\t"))["reads_subsampled"], "NA")


class SubsampleWiringTests(unittest.TestCase):
    def test_mates_use_same_seed_and_stats_count_both_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for mate in (1, 2):
                (root / f"input_{mate}.fq.gz").write_bytes(b"input")
            args = SimpleNamespace(
                sample="S1", target="50K", seed=11,
                r1=str(root / "input_1.fq.gz"), r2=str(root / "input_2.fq.gz"),
                out_r1=str(root / "S1_1.fq.gz"), out_r2=str(root / "S1_2.fq.gz"),
                stats=str(root / "S1_subsample.tsv"),
            )
            commands = []

            def sample(command, check):
                commands.append(command)
                Path(command[command.index("-o") + 1]).write_bytes(b"sampled")

            def stats(path):
                return (50_000, 500_000) if ".tmp.fq.gz" in str(path) else (100_000, 1_000_000)

            with patch.object(subsample_reads, "seqkit_stats", side_effect=stats):
                with patch.object(subsample_reads.subprocess, "run", side_effect=sample):
                    subsample_reads.run(args)
            self.assertEqual(len(commands), 2)
            for command in commands:
                self.assertEqual(command[command.index("-s") + 1], "11")
                self.assertEqual(command[command.index("-p") + 1], "0.5")
                self.assertIn("-2", command)
            with Path(args.stats).open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["reads_before"], "200000")
            self.assertEqual(row["reads_subsampled"], "100000")
            self.assertEqual(row["pairs_after"], "50000")

    def test_input_selection(self):
        source = (WORKFLOW / "Snakefile").read_text()
        block = source.split("def preprocessing_read(")[1].split("\n####", 1)[0]
        namespace = {"OUTPUT_DIR": "/out", "SANITIZE": False, "SAMPLE_TO_SUBSAMPLE": {}}
        exec("def preprocessing_read(" + block, namespace)
        read = namespace["preprocessing_read"]
        sample = SimpleNamespace(sample="S1")
        self.assertEqual(read(sample, 1), "/out/data/reads/S1_1.fq.gz")
        namespace["SAMPLE_TO_SUBSAMPLE"] = {"S1": "10Gb"}
        self.assertEqual(read(sample, 2), "/out/preprocessing/subsample/S1_2.fq.gz")
        self.assertEqual(read(SimpleNamespace(sample="S2"), 2), "/out/data/reads/S2_2.fq.gz")
        namespace["SANITIZE"] = True
        self.assertEqual(read(sample, 1), "/out/preprocessing/seqkit/S1_1.fq.gz")
        self.assertEqual(read(sample, 1, after_sanitize=False), "/out/preprocessing/subsample/S1_1.fq.gz")

    def test_profiling_read_selection(self):
        source = (WORKFLOW / "Snakefile").read_text()
        block = source.split("def profiling_read(")[1].split("def profiling_subsample_input(", 1)[0]
        namespace = {"OUTPUT_DIR": "/out", "SAMPLE_TO_SUBSAMPLE": {"S1": "50M"},
                     "PREPROCESSED_TO_READS1": {"S1": ["r1"], "S2": ["raw1"]},
                     "PREPROCESSED_TO_READS2": {"S1": ["r2"], "S2": ["raw2"]}}
        exec("def profiling_read(" + block, namespace)
        self.assertEqual(namespace["profiling_read"](SimpleNamespace(sample="S1"), 1),
                         "/out/profiling/subsample/S1_1.fq.gz")
        self.assertEqual(namespace["profiling_read"](SimpleNamespace(sample="S2"), 2), ["raw2"])
        genomes = (WORKFLOW / "rules/profiling_genomes.smk").read_text()
        self.assertEqual(genomes.count("profiling_read(wildcards, 1)"), 4)

    def test_cataloging_and_expressing_read_selection(self):
        source = (WORKFLOW / "Snakefile").read_text()
        block = source.split("def downstream_sample_reads(")[1].split("\n####", 1)[0]
        namespace = {
            "OUTPUT_DIR": "/out", "WORKFLOW": "cataloging",
            "SAMPLE_TO_SUBSAMPLE": {"S1": "10Gb"},
            "DOWNSTREAM_READS1": {"S1": ["s1_lane1", "s1_lane2"], "S2": ["s2_raw"]},
            "DOWNSTREAM_READS2": {"S1": ["s1_r2"], "S2": ["s2_r2"]},
        }
        exec("def downstream_sample_reads(" + block, namespace)
        select = namespace["downstream_sample_reads"]
        self.assertEqual(select("S1", 1), ["/out/cataloging/subsample/S1_1.fq.gz"])
        self.assertEqual(select("S2", 1), ["s2_raw"])
        namespace["WORKFLOW"] = "expressing"
        self.assertEqual(select("S1", 2), ["/out/expressing/subsample/S1_2.fq.gz"])
        catalog = (WORKFLOW / "rules/cataloging.smk").read_text()
        expressing = (WORKFLOW / "rules/expressing.smk").read_text()
        self.assertIn("downstream_sample_reads(sample, 1)", catalog)
        self.assertIn("downstream_read(wildcards, 1)", catalog)
        self.assertIn("downstream_read(wildcards, 1)", expressing)
        self.assertIn('catalog_subsample = args.subsample if args.command == "cataloging" else None',
                      (ROOT / "drakkar/cli_main.py").read_text())

    def test_both_workflow_branches_and_consumers(self):
        snakefile = (WORKFLOW / "Snakefile").read_text()
        self.assertEqual(snakefile.count('include: "rules/preprocessing_subsample.smk"'), 2)
        for name in ("preprocessing.smk", "preprocessing_ref.smk", "preprocessing_sanitize.smk"):
            source = (WORKFLOW / "rules" / name).read_text()
            self.assertIn("preprocessing_read(wildcards, 1", source)
            self.assertIn("preprocessing_read(wildcards, 2", source)
        rule = (WORKFLOW / "rules" / "preprocessing_subsample.smk").read_text()
        for field in ("input:", "output:", "log:", "threads:", "mem_mb=", "runtime=", "seqkit_module=config[\"SEQKIT_MODULE\"]"):
            self.assertIn(field, rule)


if __name__ == "__main__":
    unittest.main()
