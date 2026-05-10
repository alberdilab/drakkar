from __future__ import annotations

import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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
            self.assertIn("benchmark_jobs", metadata)
            self.assertIn("benchmark_rules", metadata)
            self.assertIn("benchmark_summary", metadata)
            self.assertEqual(
                Path(metadata["benchmark_summary"]).resolve(),
                (Path(tmpdir) / f"drakkar_{run_info['run_id']}_resources.yaml").resolve(),
            )

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

    def test_run_logging_default_shows_summary_and_usage_guide(self) -> None:
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
            self.assertIn("EXECUTION SUMMARY", output)
            self.assertIn("Status: failed", output)
            self.assertIn("HOW TO INSPECT MORE", output)
            self.assertIn("Failure excerpt or tail:", output)
            self.assertNotIn("SNAKEMAKE LOG", output)
            self.assertNotIn("Most recent failure excerpt:", output)

    def test_run_logging_excerpt_shows_failure_excerpt_from_latest_run(self) -> None:
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
                exit_code = cli_module.run_logging(tmpdir, tail=5, excerpt=True)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("SNAKEMAKE LOG", output)
            self.assertIn("Most recent failure excerpt:", output)
            self.assertIn("RuleException in rule map_reads:", output)
            self.assertIn("output: sample.bam", output)

    def test_run_logging_summary_reports_progress_rules_and_error_types(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir)
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            cli_module.update_launch_metadata(
                run_info["metadata_path"],
                status="failed",
                current_workflow="cataloging",
                exit_code=1,
            )
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "Job stats:",
                        "job            count",
                        "-----------  -------",
                        "all                1",
                        "map_reads          2",
                        "dereplicate        1",
                        "total              4",
                        "",
                        "rule map_reads:",
                        "    jobid: 1",
                        "Finished jobid: 1 (Rule: map_reads)",
                        "1 of 4 steps (25%) done",
                        "",
                        "rule map_reads:",
                        "    jobid: 2",
                        "Finished jobid: 2 (Rule: map_reads)",
                        "2 of 4 steps (50%) done",
                        "",
                        "localrule dereplicate:",
                        "    jobid: 3",
                        "RuleException in rule dereplicate:",
                        "Error in rule dereplicate:",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_logging(tmpdir, summary=True)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("EXECUTION SUMMARY", output)
            self.assertIn("Planned jobs: 4", output)
            self.assertIn("Workflow progress: 50% (2/4 steps)", output)
            self.assertIn("Rules observed: 2 unique, 3 executions", output)
            self.assertIn("Failed rules detected: 1", output)
            self.assertIn("Error types: RuleException (1), RuleError (1)", output)
            self.assertIn("HOW TO INSPECT MORE", output)
            self.assertNotIn("SNAKEMAKE LOG", output)

    def test_parse_snakemake_submitted_launches_infers_attempts_and_skips_localrules(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "rule map_reads:",
                        "    jobid: 1",
                        "    wildcards: sample=A",
                        "    threads: 8",
                        "    resources: mem_mb=16000, runtime=30, tmpdir=/tmp",
                        "Submitted job 1 with external jobid '101'.",
                        "",
                        "rule map_reads:",
                        "    jobid: 2",
                        "    wildcards: sample=A",
                        "    threads: 8",
                        "    resources: mem_mb=32000, runtime=60, tmpdir=/tmp",
                        "Submitted job 2 with external jobid '102'.",
                        "",
                        "localrule summarize:",
                        "    jobid: 3",
                        "    threads: 1",
                        "Submitted job 3 with external jobid '103'.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual(len(launches), 2)
            self.assertEqual([launch["attempt"] for launch in launches], [1, 2])
            self.assertEqual([launch["external_jobid"] for launch in launches], ["101", "102"])
            self.assertEqual(launches[0]["requested_mem_mb"], 16000.0)
            self.assertEqual(launches[1]["requested_runtime_min"], 60.0)

    def test_parse_snakemake_submitted_launches_accepts_sbatch_style_submission_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "rule map_reads:",
                        "    jobid: 1",
                        "    wildcards: sample=A",
                        "    threads: 8",
                        "    resources: mem_mb=16000, runtime=30, tmpdir=/tmp",
                        "Submitted batch job 501234",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual(len(launches), 1)
            self.assertEqual(launches[0]["internal_jobid"], "1")
            self.assertEqual(launches[0]["external_jobid"], "501234")

    def test_parse_snakemake_submitted_launches_accepts_slurm_executor_submission_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "[Sun May 10 03:27:11 2026]",
                        "rule map_reads:",
                        "    input: reads/A_R1.fastq.gz, reads/A_R2.fastq.gz",
                        "    output: mapped/A.bam",
                        "    jobid: 17",
                        "    wildcards: sample=A",
                        "    threads: 8",
                        "    resources: mem_mb=16000, runtime=30, tmpdir=<TBD>",
                        "Job 17 has been submitted with SLURM jobid 8281752 (log: .snakemake/slurm_logs/rule_map_reads/A/8281752.log).",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual(len(launches), 1)
            self.assertEqual(launches[0]["internal_jobid"], "17")
            self.assertEqual(launches[0]["external_jobid"], "8281752")
            self.assertEqual(launches[0]["rule"], "map_reads")

    def test_parse_snakemake_submitted_launches_normalizes_embedded_sbatch_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "rule map_reads:",
                        "    jobid: 1",
                        "    wildcards: sample=A",
                        "    threads: 8",
                        "    resources: mem_mb=16000, runtime=30, tmpdir=/tmp",
                        "Submitted job 1 with external jobid 'Submitted batch job 6018753'.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual(len(launches), 1)
            self.assertEqual(launches[0]["external_jobid"], "6018753")

    def test_generate_run_benchmark_writes_reports_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir, profile="slurm")
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            cli_module.update_launch_metadata(
                run_info["metadata_path"],
                status="failed",
                current_workflow="cataloging",
                exit_code=1,
            )
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "rule assemble:",
                        "    jobid: 1",
                        "    wildcards: assembly=A",
                        "    threads: 4",
                        "    resources: mem_mb=8000, runtime=30, tmpdir=/tmp",
                        "Submitted job 1 with external jobid '201'.",
                        "",
                        "rule assemble:",
                        "    jobid: 2",
                        "    wildcards: assembly=A",
                        "    threads: 4",
                        "    resources: mem_mb=16000, runtime=60, tmpdir=/tmp",
                        "Submitted job 2 with external jobid '202'.",
                        "",
                        "rule annotate:",
                        "    jobid: 3",
                        "    wildcards: mag=M1",
                        "    threads: 1",
                        "    resources: mem_mb=4000, runtime=15, tmpdir=/tmp",
                        "Submitted job 3 with external jobid '203'.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            accounting = {
                "201": {
                    "external_jobid": "201",
                    "state": "OUT_OF_MEMORY",
                    "exit_code": "0:0",
                    "elapsed_sec": 300,
                    "cpu_time_sec": 900,
                    "alloc_cpus": 4,
                    "max_rss_mb": 7900.0,
                    "timelimit_raw_min": 30,
                },
                "202": {
                    "external_jobid": "202",
                    "state": "COMPLETED",
                    "exit_code": "0:0",
                    "elapsed_sec": 600,
                    "cpu_time_sec": 2000,
                    "alloc_cpus": 4,
                    "max_rss_mb": 12000.0,
                    "timelimit_raw_min": 60,
                },
                "203": {
                    "external_jobid": "203",
                    "state": "TIMEOUT",
                    "exit_code": "0:0",
                    "elapsed_sec": 900,
                    "cpu_time_sec": 800,
                    "alloc_cpus": 1,
                    "max_rss_mb": 3500.0,
                    "timelimit_raw_min": 15,
                },
            }

            with patch.object(cli_module, "query_sacct_for_jobs", return_value=accounting):
                result = cli_module.generate_run_benchmark(tmpdir, metadata_path=run_info["metadata_path"], quiet=True)

            self.assertIsNotNone(result)
            self.assertEqual(result["status"], "generated")
            self.assertEqual(result["summary"]["benchmarked_launches"], 3)
            self.assertEqual(result["summary"]["retries"], 1)
            self.assertEqual(result["summary"]["oom_launches"], 1)
            self.assertEqual(result["summary"]["timeout_launches"], 1)

            jobs_path = Path(result["paths"]["jobs"])
            rules_path = Path(result["paths"]["rules"])
            summary_path = Path(result["paths"]["summary"])
            self.assertTrue(jobs_path.exists())
            self.assertTrue(rules_path.exists())
            self.assertTrue(summary_path.exists())
            self.assertEqual(summary_path, Path(tmpdir) / f"drakkar_{run_info['run_id']}_resources.yaml")

            jobs_text = jobs_path.read_text(encoding="utf-8")
            self.assertIn("attempt", jobs_text)
            self.assertIn("OUT_OF_MEMORY", jobs_text)
            self.assertIn("\t2\tassemble|wildcards|assembly=A\t", jobs_text)

            summary = yaml.safe_load(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "generated")
            self.assertEqual(summary["benchmarked_launches"], 3)
            self.assertEqual(summary["retries"], 1)
            self.assertIn("rules", summary)

    def test_generate_run_benchmark_writes_empty_tables_when_no_submitted_jobs_are_found(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="profiling", output=tmpdir, profile="slurm")
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "rule dereplicate:",
                        "    jobid: 1",
                        "Finished jobid: 1 (Rule: dereplicate)",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = cli_module.generate_run_benchmark(tmpdir, metadata_path=run_info["metadata_path"], quiet=True)

            self.assertEqual(result["status"], "no_submitted_jobs")
            self.assertTrue(Path(result["paths"]["jobs"]).exists())
            self.assertTrue(Path(result["paths"]["rules"]).exists())
            self.assertEqual(Path(result["paths"]["jobs"]).read_text(encoding="utf-8").strip(), "\t".join(cli_module.BENCHMARK_JOB_FIELDS))
            self.assertEqual(Path(result["paths"]["rules"]).read_text(encoding="utf-8").strip(), "\t".join(cli_module.BENCHMARK_RULE_FIELDS))

    def test_generate_run_benchmark_writes_root_status_file_for_non_slurm_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir, profile="local")
            run_info = cli_module.write_launch_metadata(args, tmpdir)

            result = cli_module.generate_run_benchmark(tmpdir, metadata_path=run_info["metadata_path"], quiet=True)

            self.assertIsNotNone(result)
            self.assertEqual(result["status"], "unsupported_profile")
            summary_path = Path(result["paths"]["summary"])
            self.assertTrue(summary_path.exists())
            self.assertEqual(summary_path, Path(tmpdir) / f"drakkar_{run_info['run_id']}_resources.yaml")

            summary = yaml.safe_load(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "unsupported_profile")
            self.assertEqual(summary["profile"], "local")

    def test_generate_run_benchmark_honors_skip_benchmark_flag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir, profile="slurm", skip_benchmark=True)
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            Path(run_info["snakemake_log_path"]).write_text("rule assemble:\n    jobid: 1\nSubmitted batch job 12345\n", encoding="utf-8")

            result = cli_module.generate_run_benchmark(tmpdir, metadata_path=run_info["metadata_path"], quiet=True)

            self.assertEqual(result["status"], "skipped")
            summary = yaml.safe_load(Path(result["paths"]["summary"]).read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "skipped")
            self.assertIn("--skip-benchmark", summary["message"])

    def test_run_logging_summary_prints_benchmark_section_for_slurm_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="cataloging", output=tmpdir, profile="slurm")
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            cli_module.update_launch_metadata(
                run_info["metadata_path"],
                status="success",
                current_workflow="cataloging",
                exit_code=0,
            )
            Path(run_info["snakemake_log_path"]).write_text(
                "\n".join(
                    [
                        "Job stats:",
                        "job            count",
                        "-----------  -------",
                        "assemble           1",
                        "total              1",
                        "",
                        "rule assemble:",
                        "    jobid: 1",
                        "    wildcards: assembly=A",
                        "    threads: 4",
                        "    resources: mem_mb=8000, runtime=30, tmpdir=/tmp",
                        "Submitted job 1 with external jobid '301'.",
                        "Finished jobid: 1 (Rule: assemble)",
                        "1 of 1 steps (100%) done",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with patch.object(
                cli_module,
                "query_sacct_for_jobs",
                return_value={
                    "301": {
                        "external_jobid": "301",
                        "state": "COMPLETED",
                        "exit_code": "0:0",
                        "elapsed_sec": 600,
                        "cpu_time_sec": 1800,
                        "alloc_cpus": 4,
                        "max_rss_mb": 6000.0,
                        "timelimit_raw_min": 30,
                    }
                },
            ):
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer):
                    exit_code = cli_module.run_logging(tmpdir, summary=True)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("RESOURCE BENCHMARK", output)
            self.assertIn("Benchmarked launches: 1", output)
            self.assertIn("Relaunches detected: 0", output)
            self.assertIn("Weighted CPU efficiency: 75.0%", output)


if __name__ == "__main__":
    unittest.main()
