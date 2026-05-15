from __future__ import annotations

import argparse
import contextlib
import io
import json
import re
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from drakkar import cli as cli_module


class StatusCommandTests(unittest.TestCase):
    def write_sample_context(self, output_dir: str) -> None:
        data_dir = Path(output_dir) / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        (data_dir / "sample_to_reads1.json").write_text(
            json.dumps({"A": ["A_R1.fq.gz"], "B": ["B_R1.fq.gz"]}),
            encoding="utf-8",
        )
        (data_dir / "assembly_to_samples.json").write_text(
            json.dumps({"A": ["A"], "B": ["B"]}),
            encoding="utf-8",
        )

    def write_preprocessing_log(self, log_path: Path) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "\n".join(
                [
                    "Building DAG of jobs...",
                    "Job stats:",
                    "job                         count",
                    "------------------------  -------",
                    "all                             1",
                    "concatenate_or_link_reads       2",
                    "fastp                           2",
                    "reference_map                   2",
                    "split_reads                     2",
                    "preprocessing_stats             1",
                    "total                          10",
                    "",
                    "rule concatenate_or_link_reads:",
                    "    jobid: 1",
                    "    wildcards: sample=A",
                    "Finished jobid: 1 (Rule: concatenate_or_link_reads)",
                    "",
                    "rule fastp:",
                    "    jobid: 2",
                    "    wildcards: sample=A",
                    "Finished jobid: 2 (Rule: fastp)",
                    "",
                    "rule fastp:",
                    "    jobid: 3",
                    "    wildcards: sample=B",
                    "",
                    "rule reference_map:",
                    "    jobid: 4",
                    "    wildcards: sample=A, reference=host",
                    "Finished jobid: 4 (Rule: reference_map)",
                    "",
                    "rule split_reads:",
                    "    jobid: 5",
                    "    wildcards: sample=A",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    def prepare_run(self, output_dir: str):
        args = argparse.Namespace(command="preprocessing", output=output_dir, profile="slurm")
        run_info = cli_module.write_launch_metadata(args, output_dir)
        cli_module.update_launch_metadata(
            run_info["metadata_path"],
            status="running",
            current_workflow="preprocessing",
        )
        self.write_sample_context(output_dir)
        self.write_preprocessing_log(Path(run_info["snakemake_log_path"]))
        return run_info

    def write_preprocessing_run(self, output_dir: str, run_id: str, status: str = "running") -> Path:
        output_path = Path(output_dir)
        metadata_path = output_path / f"drakkar_{run_id}.yaml"
        log_path = output_path / "log" / f"drakkar_{run_id}.snakemake.log"
        metadata_path.write_text(
            "\n".join(
                [
                    f"run_id: {run_id}",
                    f"timestamp: '2026-05-15T{run_id[-6:-4]}:{run_id[-4:-2]}:{run_id[-2:]}+00:00'",
                    f"started_at: '2026-05-15T{run_id[-6:-4]}:{run_id[-4:-2]}:{run_id[-2:]}+00:00'",
                    "command: preprocessing",
                    "modules:",
                    "- preprocessing",
                    f"status: {status}",
                    "current_workflow: preprocessing",
                    f"snakemake_log: {log_path.resolve()}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        self.write_preprocessing_log(log_path)
        return metadata_path

    def test_run_status_default_shows_rules_and_samples_hiding_helpers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.prepare_run(tmpdir)

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(output_dir=tmpdir)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("DRAKKAR STATUS", output)
            self.assertIn("OVERALL PROGRESS", output)
            self.assertIn("Rule jobs:", output)
            self.assertIn("RULE STATUS", output)
            self.assertIn("fastp", output)
            self.assertIn("reference_map", output)
            self.assertNotIn("concatenate_or_link_reads", output)
            self.assertIn("SAMPLE STATUS", output)
            self.assertIn("A", output)
            self.assertIn("2/3", output)
            self.assertIn("split_reads", output)

    def test_run_status_complete_rules_includes_helper_rules(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.prepare_run(tmpdir)

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(
                    output_dir=tmpdir,
                    show_complete=True,
                    view="rules",
                )

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("RULE STATUS", output)
            self.assertIn("concatenate_or_link_reads", output)
            self.assertNotIn("SAMPLE STATUS", output)

    def test_run_status_accepts_metadata_yaml_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_info = self.prepare_run(tmpdir)
            metadata_name = Path(run_info["metadata_path"]).name

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(
                    target=metadata_name,
                    output_dir=tmpdir,
                    view="rules",
                )

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn(f"Metadata file: {Path(run_info['metadata_path']).resolve()}", output)
            self.assertIn("RULE STATUS", output)
            self.assertNotIn("SAMPLE STATUS", output)

    def test_run_status_prompts_for_run_when_multiple_runs_are_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.write_sample_context(tmpdir)
            selected_metadata_path = self.write_preprocessing_run(
                tmpdir,
                "20260515-034520",
                status="failed",
            )
            self.write_preprocessing_run(tmpdir, "20260515-050000", status="running")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer), mock.patch(
                "drakkar.status.sys.stdin.isatty",
                return_value=True,
            ), mock.patch("builtins.input", return_value="2"):
                exit_code = cli_module.run_status(output_dir=tmpdir, view="rules")

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("AVAILABLE RUNS", output)
            self.assertIn("1. 20260515-050000: command=preprocessing, status=running", output)
            self.assertIn("2. 20260515-034520: command=preprocessing, status=failed", output)
            self.assertIn("Run ID: 20260515-034520", output)
            self.assertIn(f"Metadata file: {selected_metadata_path.resolve()}", output)

    def test_run_status_ignores_resource_summary_yaml_when_selecting_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_info = self.prepare_run(tmpdir)
            summary_path = Path(tmpdir) / f"drakkar_{run_info['run_id']}_resources.yaml"
            summary_path.write_text(
                "\n".join(
                    [
                        f"run_id: {run_info['run_id']}",
                        "command: preprocessing",
                        "status: no_submitted_jobs",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(output_dir=tmpdir, view="rules")

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn(f"Metadata file: {Path(run_info['metadata_path']).resolve()}", output)
            self.assertIn("Status: running", output)
            self.assertNotIn("Status: no_submitted_jobs", output)

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(
                    target=summary_path.name,
                    output_dir=tmpdir,
                    view="rules",
                )

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn(f"Metadata file: {Path(run_info['metadata_path']).resolve()}", output)

    def test_parse_snakemake_status_counts_finished_job_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "Building DAG of jobs...",
                        "Job stats:",
                        "job        count  min threads  max threads",
                        "prodigal       2            1            1",
                        "total          2            1            1",
                        "",
                        "rule prodigal:",
                        "    jobid: 1",
                        "    wildcards: mag=A_bin_1",
                        "Finished job 1.",
                        "1 of 2 steps (50%) done",
                        "",
                        "rule prodigal:",
                        "    jobid: 2",
                        "    wildcards: mag=B_bin_1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            summary = cli_module.parse_snakemake_status(log_path, metadata={"status": "running"})
            prodigal = next(rule for rule in summary["rules"] if rule["rule"] == "prodigal")
            self.assertEqual(prodigal["completed"], 1)
            self.assertEqual(prodigal["started"], 2)
            self.assertEqual(prodigal["total"], 2)
            self.assertEqual(prodigal["status"], "running")

    def test_run_status_shows_mag_status_for_annotating_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="annotating", output=tmpdir, profile="slurm")
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            cli_module.update_launch_metadata(
                run_info["metadata_path"],
                status="running",
                current_workflow="annotating",
            )

            self.write_sample_context(tmpdir)
            data_dir = Path(tmpdir) / "data"
            (data_dir / "mags_to_files.json").write_text(
                json.dumps(
                    {
                        "A_bin_1": f"{tmpdir}/cataloging/final/A/A_bin_1.fa",
                        "A_bin_2": f"{tmpdir}/cataloging/final/A/A_bin_2.fa",
                        "B_bin_1": f"{tmpdir}/cataloging/final/B/B_bin_1.fa",
                    }
                ),
                encoding="utf-8",
            )
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "Building DAG of jobs...",
                        "Job stats:",
                        "job        count",
                        "prodigal       3",
                        "total          3",
                        "",
                        "rule prodigal:",
                        "    jobid: 1",
                        "    wildcards: mag=A_bin_1",
                        "Finished job 1.",
                        "",
                        "rule prodigal:",
                        "    jobid: 2",
                        "    wildcards: mag=A_bin_2",
                        "Finished job 2.",
                        "",
                        "rule prodigal:",
                        "    jobid: 3",
                        "    wildcards: mag=B_bin_1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_status(output_dir=tmpdir, view="samples")

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("MAG STATUS", output)
            self.assertNotIn("SAMPLE STATUS", output)
            self.assertRegex(output, re.compile(r"^A_bin_1\s+done\s+1/1\s+complete$", re.MULTILINE))
            self.assertRegex(output, re.compile(r"^A_bin_2\s+done\s+1/1\s+complete$", re.MULTILINE))
            self.assertRegex(output, re.compile(r"^B_bin_1\s+running\s+0/1\s+prodigal$", re.MULTILINE))


if __name__ == "__main__":
    unittest.main()
