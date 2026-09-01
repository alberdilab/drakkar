from __future__ import annotations

import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
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
                Path(tmpdir) / "logging" / f"drakkar_{run_info['run_id']}.snakemake.log",
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
                (
                    Path(tmpdir) / "logging" / "benchmark" / f"drakkar_{run_info['run_id']}.resources.yaml"
                ).resolve(),
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

    def test_parse_snakemake_submitted_launches_recovers_rules_that_print_only_a_message(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "[Sun May 10 03:27:11 2026]",
                        "Job 42: Assembling EHA05984...",
                        "Job 42 has been submitted with SLURM jobid 44448391 "
                        "(log: /work/.snakemake/slurm_logs/rule_assembly/EHA05984/44448391.log).",
                        "[Sun May 10 04:11:02 2026]",
                        "Job 51: Mapping S01 reads to assembly EHA05984...",
                        "Job 51 has been submitted with SLURM jobid 44448392 "
                        "(log: /work/.snakemake/slurm_logs/rule_assembly_map/EHA05984/S01/44448392.log).",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual([launch["rule"] for launch in launches], ["assembly", "assembly_map"])
            self.assertEqual([launch["wildcards"] for launch in launches], ["EHA05984", "EHA05984,S01"])
            self.assertEqual(
                [launch["external_jobid"] for launch in launches], ["44448391", "44448392"]
            )
            self.assertIsNone(launches[0]["requested_mem_mb"])

    def test_parse_snakemake_submitted_launches_counts_retries_of_the_same_internal_job(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "drakkar_test.snakemake.log"
            log_path.write_text(
                "\n".join(
                    [
                        "Job 42: Binning contigs from assembly A using metabat2...",
                        "Job 42 has been submitted with SLURM jobid 901 "
                        "(log: .snakemake/slurm_logs/rule_metabat2/A/901.log).",
                        "Job 42: Binning contigs from assembly A using metabat2...",
                        "Job 42 has been submitted with SLURM jobid 902 "
                        "(log: .snakemake/slurm_logs/rule_metabat2/A/902.log).",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            launches = cli_module.parse_snakemake_submitted_launches(log_path)

            self.assertEqual(len(launches), 2)
            self.assertEqual([launch["attempt"] for launch in launches], [1, 2])
            self.assertEqual([launch["rule"] for launch in launches], ["metabat2", "metabat2"])

    def test_benchmark_job_row_falls_back_to_accounting_for_requested_resources(self) -> None:
        launch = {
            "launch_index": 1,
            "rule": "assembly",
            "attempt": 1,
            "logical_job_key": "assembly|wildcards|A",
            "internal_jobid": "42",
            "external_jobid": "901",
            "wildcards": "A",
            "requested_cpus": None,
            "requested_mem_mb": None,
            "requested_runtime_min": None,
        }
        accounting = {
            "state": "COMPLETED",
            "exit_code": "0:0",
            "alloc_cpus": 8,
            "elapsed_sec": 600,
            "cpu_time_sec": 4800,
            "max_rss_mb": 4000.0,
            "timelimit_raw_min": 30,
            "req_cpus": 8,
            "req_mem_mb": 8000.0,
        }

        row = cli_module.benchmark_job_row(launch, accounting)

        self.assertEqual(row["requested_cpus"], 8)
        self.assertEqual(row["requested_mem_mb"], 8000.0)
        self.assertEqual(row["requested_runtime_min"], 30)
        self.assertAlmostEqual(row["memory_efficiency"], 0.5)
        self.assertAlmostEqual(row["runtime_efficiency"], 600 / 1800)

    def test_parse_sacct_output_reads_requested_cpus_and_memory(self) -> None:
        rows = cli_module._benchmark._parse_sacct_output(
            "901|COMPLETED|0:0|600|01:20:00|8|4000M|30|8|8000M\n"
            "902|COMPLETED|0:0|600|01:20:00|8|4000M|30|8|1000Mc\n"
        )

        self.assertEqual(rows["901"]["req_cpus"], 8)
        self.assertEqual(rows["901"]["req_mem_mb"], 8000.0)
        # A per-CPU request only becomes the job total once scaled by its CPUs.
        self.assertEqual(rows["902"]["req_mem_mb"], 8000.0)

    def test_parse_sacct_output_reads_consumed_not_reserved_cpu_time(self) -> None:
        # sacct is asked for TotalCPU, which arrives as a duration rather than
        # as raw seconds. CPUTimeRAW would have been AllocCPUS x Elapsed here
        # — 4,800 s — and would report the job as 100% CPU efficient whatever
        # it actually did.
        rows = cli_module._benchmark._parse_sacct_output(
            "901|COMPLETED|0:0|600|01:20:00|8|4000M|30|8|8000M\n"
        )

        self.assertEqual(rows["901"]["cpu_time_sec"], 4800.0)
        self.assertEqual(rows["901"]["elapsed_sec"], 600)
        self.assertEqual(rows["901"]["alloc_cpus"], 8)

    def test_sacct_is_queried_for_total_cpu(self) -> None:
        with patch.object(cli_module._benchmark.subprocess, "run") as run:
            run.return_value = SimpleNamespace(returncode=0, stdout="")
            cli_module._benchmark.query_sacct_for_jobs(["901"])

        command = run.call_args[0][0]
        self.assertIn("TotalCPU", command[-1])
        self.assertNotIn("CPUTimeRAW", command[-1])

    def test_parse_slurm_duration_covers_the_formats_sacct_writes(self) -> None:
        parse = cli_module._benchmark.parse_slurm_duration_to_seconds
        # Two colon-separated fields are minutes and seconds, not hours.
        self.assertAlmostEqual(parse("12:30.500"), 750.5)
        self.assertEqual(parse("01:20:00"), 4800)
        self.assertEqual(parse("2-03:04:05"), 2 * 86400 + 3 * 3600 + 4 * 60 + 5)
        self.assertEqual(parse("00:00:00"), 0)
        for blank in ("", None, "Unknown", "not a duration"):
            self.assertIsNone(parse(blank))

    def test_cpu_efficiency_is_the_share_of_the_reservation_actually_burned(self) -> None:
        launch = {
            "launch_index": 1,
            "rule": "assembly",
            "attempt": 1,
            "logical_job_key": "assembly|wildcards|A",
            "internal_jobid": "42",
            "external_jobid": "901",
            "wildcards": "A",
            "requested_cpus": 8,
            "requested_mem_mb": 8000.0,
            "requested_runtime_min": 30,
        }
        # 8 CPUs held for 600 s is a 4,800 CPU-second reservation; the job
        # burned 1,200 of them, which is a quarter of what it held.
        accounting = {
            "state": "COMPLETED",
            "exit_code": "0:0",
            "alloc_cpus": 8,
            "elapsed_sec": 600,
            "cpu_time_sec": 1200.0,
            "max_rss_mb": 4000.0,
            "timelimit_raw_min": 30,
            "req_cpus": 8,
            "req_mem_mb": 8000.0,
        }

        row = cli_module.benchmark_job_row(launch, accounting)

        self.assertAlmostEqual(row["cpu_efficiency"], 0.25)

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
            self.assertEqual(
                summary_path,
                Path(tmpdir) / "logging" / "benchmark" / f"drakkar_{run_info['run_id']}.resources.yaml",
            )

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
            self.assertEqual(
                summary_path,
                Path(tmpdir) / "logging" / "benchmark" / f"drakkar_{run_info['run_id']}.resources.yaml",
            )

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


class LegacyLayoutTests(unittest.TestCase):
    """Output directories written before the `logging/` directory existed.

    Their run metadata sits in the output root, their Snakemake log in `log/`,
    and their benchmark roll-up and failure table use an underscore before
    `resources`/`failures`. Those directories are still read, and anything
    regenerated for such a run is written beside the files it already has
    rather than split across both layouts.
    """

    def legacy_run(self, root: Path, run_id: str = "20260503-101530") -> Path:
        """Write a legacy-layout run and return its metadata path."""
        metadata_path = root / f"drakkar_{run_id}.yaml"
        metadata_path.write_text(
            yaml.safe_dump(
                {
                    "run_id": run_id,
                    "command": "cataloging",
                    "status": "success",
                    "profile": "local",
                },
                sort_keys=False,
            ),
            encoding="utf-8",
        )
        log_path = root / "log" / f"drakkar_{run_id}.snakemake.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("Building DAG of jobs...\n", encoding="utf-8")
        return metadata_path

    def test_discover_run_metadata_finds_runs_in_the_output_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata_path = self.legacy_run(root)

            runs = cli_module.discover_run_metadata(root)

            self.assertEqual([path for path, _ in runs], [metadata_path])

    def test_discover_run_metadata_prefers_the_logging_directory(self) -> None:
        # The same run id in both layouts is one run, reported once, from the
        # current location.
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self.legacy_run(root)
            current = root / "logging" / "drakkar_20260503-101530.yaml"
            current.parent.mkdir(parents=True, exist_ok=True)
            current.write_text(
                yaml.safe_dump({"run_id": "20260503-101530", "command": "cataloging"}),
                encoding="utf-8",
            )

            runs = cli_module.discover_run_metadata(root)

            self.assertEqual([path for path, _ in runs], [current])

    def test_discover_run_metadata_reports_both_layouts_newest_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            older = self.legacy_run(root, run_id="20260503-101530")
            newer = root / "logging" / "drakkar_20260829-034907.yaml"
            newer.parent.mkdir(parents=True, exist_ok=True)
            newer.write_text(
                yaml.safe_dump({"run_id": "20260829-034907", "command": "profiling"}),
                encoding="utf-8",
            )

            runs = cli_module.discover_run_metadata(root)

            self.assertEqual([path for path, _ in runs], [newer, older])

    def test_find_snakemake_log_falls_back_to_the_log_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self.legacy_run(root)

            found = cli_module.find_snakemake_log(root, "20260503-101530")

            self.assertEqual(found, root / "log" / "drakkar_20260503-101530.snakemake.log")

    def test_benchmark_for_a_legacy_run_stays_in_the_legacy_layout(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata_path = self.legacy_run(root)

            result = cli_module.generate_run_benchmark(root, metadata_path=metadata_path, quiet=True)

            self.assertEqual(result["status"], "unsupported_profile")
            summary_path = Path(result["paths"]["summary"])
            self.assertEqual(summary_path, root / "drakkar_20260503-101530_resources.yaml")
            self.assertTrue(summary_path.exists())
            self.assertFalse((root / "logging").exists())

    def test_failure_report_for_a_legacy_run_stays_in_the_legacy_layout(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata_path = self.legacy_run(root)
            (root / "log" / "drakkar_20260503-101530.snakemake.log").write_text(
                "\n".join(
                    [
                        "[Mon May  3 10:15:30 2026]",
                        "rule assemble:",
                        "    jobid: 1",
                        "    wildcards: assembly=A",
                        "Error in rule assemble:",
                        "    jobid: 1",
                        "    output: cataloging/assembly/A.fna",
                        "Exiting because a job execution failed.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            report = cli_module.generate_failure_report(root, metadata_path=metadata_path)

            self.assertIsNotNone(report)
            self.assertEqual(
                Path(report["report_path"]), root / "drakkar_20260503-101530_failures.tsv"
            )
            self.assertTrue(Path(report["report_path"]).exists())
            self.assertFalse((root / "logging").exists())

    def test_job_logs_are_found_in_both_layouts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            self.legacy_run(root)
            legacy_job_log = root / "log" / "annotating" / "prodigal" / "mag1.log"
            legacy_job_log.parent.mkdir(parents=True, exist_ok=True)
            legacy_job_log.write_text("legacy\n", encoding="utf-8")
            current_job_log = root / "logging" / "annotating" / "prodigal" / "mag2.log"
            current_job_log.parent.mkdir(parents=True, exist_ok=True)
            current_job_log.write_text("current\n", encoding="utf-8")

            found = cli_module.discover_job_logs(root)

            self.assertIn(legacy_job_log, found)
            self.assertIn(current_job_log, found)

    def test_run_artefacts_are_not_reported_as_job_logs(self) -> None:
        # The metadata, log, failure table and benchmark files sit at the top of
        # the logging directory and under `benchmark/`; only module
        # subdirectories hold per-rule job logs.
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            args = argparse.Namespace(command="cataloging", output=str(root))
            run_info = cli_module.write_launch_metadata(args, str(root))
            Path(run_info["snakemake_log_path"]).write_text("log\n", encoding="utf-8")
            benchmark_table = root / "logging" / "benchmark" / "drakkar_x.jobs.tsv"
            benchmark_table.parent.mkdir(parents=True, exist_ok=True)
            benchmark_table.write_text("rule\n", encoding="utf-8")

            self.assertEqual(cli_module.discover_job_logs(root), [])


if __name__ == "__main__":
    unittest.main()
