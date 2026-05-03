from __future__ import annotations

import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import yaml

from drakkar import cli as cli_module


class LoggingCommandTests(unittest.TestCase):
    def test_write_launch_metadata_records_run_id_and_snakemake_log_for_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir)

            run_info = cli_module.write_launch_metadata(args, tmpdir)

            self.assertIsNotNone(run_info)
            metadata_path = Path(run_info["metadata_path"])
            self.assertTrue(metadata_path.exists())
            self.assertEqual(metadata_path.name, f"drakkar_{run_info['run_id']}.yaml")
            self.assertEqual(
                Path(run_info["snakemake_log_path"]),
                Path(tmpdir) / "log" / f"drakkar_{run_info['run_id']}.snakemake.log",
            )

            metadata = yaml.safe_load(metadata_path.read_text(encoding="utf-8"))
            self.assertEqual(metadata["run_id"], run_info["run_id"])
            self.assertEqual(metadata["status"], "prepared")
            self.assertEqual(metadata["command"], "cataloging")
            self.assertIn("snakemake_log", metadata)

    def test_run_subprocess_with_logging_updates_metadata_and_writes_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir)
            run_info = cli_module.write_launch_metadata(args, tmpdir)

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_module.run_subprocess_with_logging(
                    ["/bin/sh", "-c", "printf 'alpha\\nbeta\\n'"],
                    run_info=run_info,
                    workflow_name="cataloging",
                )

            log_path = Path(run_info["snakemake_log_path"])
            self.assertEqual(log_path.read_text(encoding="utf-8"), "alpha\nbeta\n")
            metadata = yaml.safe_load(Path(run_info["metadata_path"]).read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "success")
            self.assertEqual(metadata["exit_code"], 0)
            self.assertEqual(metadata["current_workflow"], "cataloging")
            self.assertIn("alpha", buffer.getvalue())
            self.assertIn("beta", buffer.getvalue())

    def test_run_logging_shows_failure_excerpt_from_latest_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="profiling", output=tmpdir)
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            cli_module.update_launch_metadata(
                run_info["metadata_path"],
                status="failed",
                current_workflow="profiling",
                exit_code=1,
            )
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "Building DAG of jobs...",
                        "Something happened",
                        "RuleException in rule map_reads:",
                        "jobid: 7",
                        "output: sample.bam",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_logging(tmpdir, tail=5)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("RUN SUMMARY", output)
            self.assertIn("Status: failed", output)
            self.assertIn("Most recent failure excerpt:", output)
            self.assertIn("RuleException in rule map_reads:", output)
            self.assertIn("output: sample.bam", output)


if __name__ == "__main__":
    unittest.main()
