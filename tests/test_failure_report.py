from __future__ import annotations

import argparse
import contextlib
import csv
import io
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml

from drakkar import cli as cli_module
from drakkar import failures as failures_module


SLURM_FAILURE_LOG = """Building DAG of jobs...
Job stats:
job                count
---------------  -------
samtools_stats         1
singlem                2
total                  3

[Fri Aug 14 01:10:00 2026]
rule samtools_stats:
    input: preprocessing/bowtie2/PR04534.bam
    output: preprocessing/samtools/PR04534.stats.txt
    jobid: 55
    wildcards: sample=PR04534
    resources: mem_mb=8000, runtime=60

Job 55 has been submitted with SLURM jobid 45360507 (log: JOBLOG).

[Fri Aug 14 01:22:12 2026]
Error in rule samtools_stats:
    message: SLURM-job '45360507' failed, SLURM status is: 'FAILED'. For further error details see the cluster/cloud log and the log files of the involved rule(s).
    jobid: 55
    input: preprocessing/bowtie2/PR04534.bam
    output: preprocessing/samtools/PR04534.stats.txt
    log: JOBLOG (check log file(s) for error details)
    shell:

        samtools index preprocessing/bowtie2/PR04534.bam

        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    external_jobid: 45360507

Trying to restart job 41.

[Fri Aug 14 01:10:05 2026]
rule singlem:
    input: preprocessing/final/PR04533_1.fq.gz
    output: preprocessing/singlem/PR04533_cond.tsv
    jobid: 41
    wildcards: sample=PR04533
    resources: mem_mb=36696, runtime=72

[Fri Aug 14 01:25:03 2026]
Error in rule singlem:
    message: SLURM-job '45359940' failed, SLURM status is: 'TIMEOUT'. For further error details see the cluster/cloud log and the log files of the involved rule(s).
    jobid: 41
    output: preprocessing/singlem/PR04533_cond.tsv
    external_jobid: 45359940

[Fri Aug 14 01:30:03 2026]
rule singlem:
    output: preprocessing/singlem/PR04535_cond.tsv
    jobid: 37
    wildcards: sample=PR04535
    resources: mem_mb=36696, runtime=72

[Fri Aug 14 01:40:03 2026]
Error in rule singlem:
    message: SLURM-job '45360510' failed, SLURM status is: 'OUT_OF_MEMORY'. For further error details see the cluster/cloud log and the log files of the involved rule(s).
    jobid: 37
    output: preprocessing/singlem/PR04535_cond.tsv
    external_jobid: 45360510

Finished job 37.
Exiting because a job execution failed. Look above for error message
"""

LOCAL_FAILURE_LOG = """Building DAG of jobs...
[Fri Aug 14 02:00:00 2026]
rule fastp:
    input: reads/A_1.fq.gz
    output: preprocessing/fastp/A_1.fq.gz
    jobid: 3
    wildcards: sample=A

ERROR: sequence file A_1.fq.gz is malformed: unexpected end of file
[Fri Aug 14 02:01:00 2026]
Error in rule fastp:
    jobid: 3
    input: reads/A_1.fq.gz
    output: preprocessing/fastp/A_1.fq.gz
    shell:

        fastp -i reads/A_1.fq.gz -o preprocessing/fastp/A_1.fq.gz
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Building DAG of jobs...
MissingInputException in rule bowtie2 in file /pkg/workflow/preprocessing.smk, line 42:
Missing input files for rule bowtie2:
    affected files:
        reference/host.fna.gz

"""


class FailureReportTests(unittest.TestCase):
    def write_slurm_run(self, tmpdir: str) -> tuple[Path, Path]:
        output_path = Path(tmpdir)
        job_log = output_path / ".snakemake" / "slurm_logs" / "rule_samtools_stats" / "45360507.log"
        job_log.parent.mkdir(parents=True, exist_ok=True)
        job_log.write_text(
            "\n".join(
                [
                    "[E::hts_open_format] Failed to open file PR04534.bam : No such file or directory",
                    'samtools index: failed to open "PR04534.bam": No such file or directory',
                ]
            ),
            encoding="utf-8",
        )
        log_path = output_path / "log" / "drakkar_20260814-011000.snakemake.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(SLURM_FAILURE_LOG.replace("JOBLOG", str(job_log)), encoding="utf-8")
        metadata_path = output_path / "drakkar_20260814-011000.yaml"
        metadata_path.write_text(
            yaml.safe_dump(
                {
                    "run_id": "20260814-011000",
                    "command": "preprocessing",
                    "status": "failed",
                    "snakemake_log": str(log_path),
                    "arguments": {"profile": "slurm"},
                }
            ),
            encoding="utf-8",
        )
        return log_path, metadata_path

    def test_parse_snakemake_failures_extracts_rule_jobid_and_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path, _ = self.write_slurm_run(tmpdir)

            parsed = failures_module.parse_snakemake_failures(log_path)

            self.assertTrue(parsed["has_log"])
            self.assertEqual(len(parsed["failures"]), 3)
            first = parsed["failures"][0]
            self.assertEqual(first["rule"], "samtools_stats")
            self.assertEqual(first["internal_jobid"], "55")
            self.assertEqual(first["external_jobid"], "45360507")
            self.assertEqual(first["slurm_state"], "FAILED")
            self.assertEqual(first["wildcards"], "sample=PR04534")
            self.assertFalse(first["recovered"])
            self.assertTrue(parsed["failures"][2]["recovered"])

    def test_build_failure_report_classifies_and_recommends_actions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path, _ = self.write_slurm_run(tmpdir)

            report = failures_module.build_failure_report(log_path, run_id="20260814-011000")

            by_target = {row["target"]: row for row in report["rows"]}
            self.assertEqual(by_target["PR04534"]["category"], "command-error")
            self.assertIn("samtools index", by_target["PR04534"]["detail"])
            self.assertEqual(by_target["PR04533"]["category"], "timeout")
            self.assertIn("72 min", by_target["PR04533"]["detail"])
            self.assertIn("--time-multiplier", by_target["PR04533"]["action"])
            self.assertEqual(by_target["PR04535"]["category"], "out-of-memory")
            self.assertEqual(by_target["PR04535"]["status"], "recovered")
            self.assertEqual(report["summary"]["failed_jobs"], 2)
            self.assertEqual(report["summary"]["failed_rules"], 2)
            self.assertEqual(report["summary"]["recovered_jobs"], 1)
            self.assertEqual(report["summary"]["verdict"], "needs-changes")

    def test_resource_failures_only_yield_rerun_with_more_resources_verdict(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "run.log"
            log_path.write_text(
                "\n".join(
                    [
                        "Building DAG of jobs...",
                        "rule singlem:",
                        "    jobid: 41",
                        "    wildcards: sample=PR04533",
                        "    resources: mem_mb=36696, runtime=72",
                        "",
                        "Error in rule singlem:",
                        "    message: SLURM-job '1' failed, SLURM status is: 'TIMEOUT'.",
                        "    jobid: 41",
                        "    external_jobid: 1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            report = failures_module.build_failure_report(log_path)

            self.assertEqual(report["summary"]["verdict"], "rerun-with-more-resources")

    def test_local_failures_use_preceding_stderr_and_report_workflow_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "run.log"
            log_path.write_text(LOCAL_FAILURE_LOG, encoding="utf-8")

            report = failures_module.build_failure_report(log_path)

            self.assertEqual(len(report["rows"]), 1)
            row = report["rows"][0]
            self.assertEqual(row["rule"], "fastp")
            self.assertEqual(row["target"], "A")
            self.assertEqual(row["category"], "command-error")
            self.assertIn("malformed", row["detail"])
            self.assertEqual(len(report["workflow_errors"]), 1)
            self.assertEqual(report["workflow_errors"][0]["category"], "missing-input")

    def test_generate_failure_report_writes_table_and_updates_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            _, metadata_path = self.write_slurm_run(tmpdir)

            report = failures_module.generate_failure_report(tmpdir, metadata_path=metadata_path)

            report_path = Path(report["report_path"])
            self.assertTrue(report_path.exists())
            self.assertEqual(report_path.name, "drakkar_20260814-011000.failures.tsv")
            with open(report_path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["rule"], "samtools_stats")
            self.assertEqual(rows[0]["external_jobid"], "45360507")
            metadata = yaml.safe_load(Path(metadata_path).read_text(encoding="utf-8"))
            self.assertEqual(metadata["failed_jobs"], 2)
            self.assertEqual(metadata["failure_verdict"], "needs-changes")
            self.assertEqual(Path(metadata["failure_report"]), report_path)

    def test_logging_command_prints_failure_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.write_slurm_run(tmpdir)

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = cli_module.run_logging(tmpdir, failures=True)

            output = buffer.getvalue()
            self.assertEqual(exit_code, 0)
            self.assertIn("FAILURE REPORT", output)
            self.assertIn("samtools_stats", output)
            self.assertIn("timeout", output)
            self.assertIn("WHAT TO DO NEXT", output)
            self.assertIn("--time-multiplier", output)
            self.assertIn("FAILED JOB DETAILS", output)

    def test_failed_workflow_run_prints_failure_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(command="preprocessing", output=tmpdir, profile="slurm")
            run_info = cli_module.write_launch_metadata(args, tmpdir)
            script = (
                "printf 'Building DAG of jobs...\\n"
                "rule fastp:\\n    jobid: 3\\n    wildcards: sample=A\\n\\n"
                "Error in rule fastp:\\n    jobid: 3\\n    output: fastp/A.fq.gz\\n\\n'; exit 1"
            )

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                with self.assertRaises(subprocess.CalledProcessError):
                    cli_module.run_subprocess_with_logging(
                        ["/bin/sh", "-c", script],
                        run_info=run_info,
                        workflow_name="preprocessing",
                    )

            output = buffer.getvalue()
            self.assertIn("FAILURE REPORT", output)
            self.assertIn("fastp", output)
            self.assertIn("Verdict:", output)
            self.assertTrue((Path(tmpdir) / "log" / f"drakkar_{run_info['run_id']}.failures.tsv").exists())


if __name__ == "__main__":
    unittest.main()
