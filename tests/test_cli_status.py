from __future__ import annotations

import argparse
import contextlib
import io
import json
import tempfile
import unittest
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


if __name__ == "__main__":
    unittest.main()
