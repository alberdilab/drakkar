import csv
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import yaml

from drakkar.cli_context import INFO, RESET
from drakkar.output import print
from drakkar.run_metadata import build_benchmark_paths, build_snakemake_log_path, load_metadata_file, update_launch_metadata

def get_run_profile(metadata):
    arguments = (metadata or {}).get("arguments") or {}
    return str(arguments.get("profile", "")).strip().lower() or None

def should_benchmark_run(metadata):
    return get_run_profile(metadata) == "slurm"

def should_skip_benchmark(metadata):
    arguments = {}
    if isinstance(metadata, dict):
        arguments = metadata.get("arguments") or {}
    return bool(arguments.get("skip_benchmark"))

def parse_int_or_none(value):
    if value in (None, ""):
        return None
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None

def parse_float_or_none(value):
    if value in (None, ""):
        return None
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None

def parse_resource_assignments(text):
    resources = {}
    if not text:
        return resources
    for item in str(text).split(","):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        resources[key.strip()] = value.strip()
    return resources

def parse_slurm_job_id(value):
    text = str(value or "").strip().strip("'\"").rstrip(".")
    if not text:
        return None

    simple_match = re.fullmatch(r"\d+(?:_\d+)?(?:\.\d+)?(?:;[A-Za-z0-9_.-]+)?", text)
    if simple_match:
        return text.split(";", 1)[0]

    matches = re.findall(r"(?<![\w.])\d+(?:_\d+)?(?:\.\d+)?(?:;[A-Za-z0-9_.-]+)?(?![\w.])", text)
    if matches:
        return matches[-1].split(";", 1)[0]
    return text

def parse_slurm_memory_to_mb(value):
    text = str(value or "").strip()
    if not text or text in {"Unknown", "None", "N/A"}:
        return None
    if text[-1:].lower() in {"c", "n"}:
        text = text[:-1]
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([kmgtpe]?)", text, flags=re.IGNORECASE)
    if not match:
        return None
    amount = float(match.group(1))
    unit = match.group(2).upper()
    factors = {
        "": 1,
        "K": 1 / 1024,
        "M": 1,
        "G": 1024,
        "T": 1024 * 1024,
        "P": 1024 * 1024 * 1024,
        "E": 1024 * 1024 * 1024 * 1024,
    }
    factor = factors.get(unit)
    if factor is None:
        return None
    return round(amount * factor, 3)

def safe_ratio(numerator, denominator):
    if numerator in (None, "") or denominator in (None, "", 0):
        return None
    try:
        denominator_value = float(denominator)
        if denominator_value == 0:
            return None
        return float(numerator) / denominator_value
    except (TypeError, ValueError, ZeroDivisionError):
        return None

def median_or_none(values):
    filtered = [float(value) for value in values if value not in (None, "")]
    if not filtered:
        return None
    return statistics.median(filtered)

def format_hours(value):
    if value is None:
        return "unknown"
    return f"{value:.2f}"

def format_megabytes(value):
    if value is None:
        return "unknown"
    if value >= 1024:
        return f"{value / 1024:.2f} GB"
    return f"{value:.0f} MB"

def format_percent(value):
    if value is None:
        return "unknown"
    return f"{value * 100:.1f}%"

def build_logical_job_key(launch):
    for field in ("wildcards", "output", "input"):
        value = str(launch.get(field, "")).strip()
        if value:
            return f"{launch.get('rule', 'unknown')}|{field}|{value}"
    return f"{launch.get('rule', 'unknown')}|singleton"

def parse_snakemake_submitted_launches(log_path):
    path = Path(log_path)
    if not path.exists():
        return []

    block_by_jobid = {}
    current_block = None
    launches = []
    launch_order = 0
    launched_internal_jobids = set()

    def register_launch(block, internal_jobid, external_jobid):
        nonlocal launch_order
        external_jobid = parse_slurm_job_id(external_jobid)
        if not block or block.get("local") or not internal_jobid or not external_jobid:
            return
        if internal_jobid in launched_internal_jobids:
            return

        launch_order += 1
        resources = block.get("resources") or {}
        requested_mem_mb = parse_float_or_none(resources.get("mem_mb"))
        if requested_mem_mb is None:
            requested_mem_mb = parse_float_or_none(resources.get("mem_mib"))
        requested_runtime_min = parse_float_or_none(resources.get("runtime"))
        launch = {
            "launch_index": launch_order,
            "rule": block.get("rule"),
            "rule_type": block.get("rule_type"),
            "internal_jobid": internal_jobid,
            "external_jobid": str(external_jobid).rstrip("."),
            "threads": block.get("threads"),
            "requested_cpus": block.get("threads"),
            "requested_mem_mb": requested_mem_mb,
            "requested_runtime_min": requested_runtime_min,
            "wildcards": block.get("wildcards", ""),
            "input": block.get("input", ""),
            "output": block.get("output", ""),
            "resources": resources,
        }
        launch["logical_job_key"] = build_logical_job_key(launch)
        launches.append(launch)
        launched_internal_jobids.add(internal_jobid)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue

            block_match = re.match(r"^(localrule|localcheckpoint|rule|checkpoint)\s+(.+?):\s*$", stripped)
            if block_match:
                current_block = {
                    "rule_type": block_match.group(1),
                    "rule": block_match.group(2),
                    "local": block_match.group(1).startswith("local"),
                    "internal_jobid": None,
                    "threads": None,
                    "wildcards": "",
                    "input": "",
                    "output": "",
                    "resources_raw": "",
                    "resources": {},
                }
                continue

            if current_block is not None:
                match = re.match(r"^jobid:\s*(\d+)\s*$", stripped)
                if match:
                    current_block["internal_jobid"] = match.group(1)
                    block_by_jobid[match.group(1)] = current_block
                    continue
                if stripped.startswith("threads:"):
                    current_block["threads"] = parse_int_or_none(stripped.split(":", 1)[1])
                    continue
                if stripped.startswith("wildcards:"):
                    current_block["wildcards"] = stripped.split(":", 1)[1].strip()
                    continue
                if stripped.startswith("input:"):
                    current_block["input"] = stripped.split(":", 1)[1].strip()
                    continue
                if stripped.startswith("output:"):
                    current_block["output"] = stripped.split(":", 1)[1].strip()
                    continue
                if stripped.startswith("resources:"):
                    current_block["resources_raw"] = stripped.split(":", 1)[1].strip()
                    current_block["resources"] = parse_resource_assignments(current_block["resources_raw"])
                    continue

            submission_match = re.search(
                r"Submitted (?:group )?job(?:id)?\s+(\d+).*?external jobid\s+(.+?)\.?$",
                stripped,
            )
            if submission_match:
                internal_jobid = submission_match.group(1)
                external_jobid = submission_match.group(2).strip().strip("'\"")
                register_launch(
                    block_by_jobid.get(internal_jobid),
                    internal_jobid,
                    external_jobid,
                )
                continue

            slurm_executor_submission_match = re.search(
                r"^Job\s+(\S+)\s+has been submitted with SLURM jobid\s+(.+?)(?:\s+\(log:|$)",
                stripped,
            )
            if slurm_executor_submission_match:
                internal_jobid = slurm_executor_submission_match.group(1)
                block = block_by_jobid.get(internal_jobid)
                if (
                    block is None
                    and current_block is not None
                    and current_block.get("internal_jobid") == internal_jobid
                ):
                    block = current_block
                register_launch(
                    block,
                    internal_jobid,
                    slurm_executor_submission_match.group(2),
                )
                continue

            batch_submission_match = re.search(
                r"Submitted batch job\s+['\"]?([^'\"\s]+)['\"]?\.?$",
                stripped,
            )
            if batch_submission_match and current_block is not None:
                register_launch(
                    current_block,
                    current_block.get("internal_jobid"),
                    batch_submission_match.group(1),
                )
                continue

            external_only_match = re.search(
                r"Submitted .*?external jobid\s+(.+?)\.?$",
                stripped,
            )
            if external_only_match and current_block is not None:
                external_jobid = external_only_match.group(1).strip().strip("'\"")
                register_launch(
                    current_block,
                    current_block.get("internal_jobid"),
                    external_jobid,
                )

    attempt_counter = Counter()
    for launch in launches:
        attempt_counter[launch["logical_job_key"]] += 1
        launch["attempt"] = attempt_counter[launch["logical_job_key"]]

    return launches

def query_sacct_for_jobs(job_ids):
    if not job_ids:
        return {}

    command = [
        "sacct",
        "-X",
        "-P",
        "-n",
        "--units=M",
        "-j",
        ",".join(str(job_id) for job_id in job_ids),
        "--format",
        "JobIDRaw,State,ExitCode,ElapsedRaw,CPUTimeRAW,AllocCPUS,MaxRSS,TimelimitRaw",
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None

    rows = {}
    for raw_line in result.stdout.splitlines():
        if not raw_line.strip():
            continue
        parts = raw_line.split("|")
        if len(parts) != 8:
            continue
        job_id = parts[0].strip()
        if not job_id:
            continue
        rows[job_id] = {
            "external_jobid": job_id,
            "state": parts[1].strip(),
            "exit_code": parts[2].strip(),
            "elapsed_sec": parse_int_or_none(parts[3]),
            "cpu_time_sec": parse_int_or_none(parts[4]),
            "alloc_cpus": parse_int_or_none(parts[5]),
            "max_rss_mb": parse_slurm_memory_to_mb(parts[6]),
            "timelimit_raw_min": parse_int_or_none(parts[7]),
        }
    return rows

def benchmark_job_row(launch, accounting_row):
    row = {
        "launch_index": launch.get("launch_index"),
        "rule": launch.get("rule"),
        "attempt": launch.get("attempt"),
        "logical_job_key": launch.get("logical_job_key"),
        "internal_jobid": launch.get("internal_jobid"),
        "external_jobid": launch.get("external_jobid"),
        "wildcards": launch.get("wildcards", ""),
        "requested_cpus": launch.get("requested_cpus"),
        "requested_mem_mb": launch.get("requested_mem_mb"),
        "requested_runtime_min": launch.get("requested_runtime_min"),
        "state": None,
        "exit_code": None,
        "alloc_cpus": None,
        "elapsed_sec": None,
        "cpu_time_sec": None,
        "max_rss_mb": None,
        "cpu_efficiency": None,
        "memory_efficiency": None,
        "runtime_efficiency": None,
        "oom": False,
        "timeout": False,
    }
    if accounting_row:
        row.update(
            {
                "state": accounting_row.get("state"),
                "exit_code": accounting_row.get("exit_code"),
                "alloc_cpus": accounting_row.get("alloc_cpus"),
                "elapsed_sec": accounting_row.get("elapsed_sec"),
                "cpu_time_sec": accounting_row.get("cpu_time_sec"),
                "max_rss_mb": accounting_row.get("max_rss_mb"),
            }
        )
        alloc_cpus = row.get("alloc_cpus")
        elapsed_sec = row.get("elapsed_sec")
        cpu_time_sec = row.get("cpu_time_sec")
        if alloc_cpus and elapsed_sec:
            row["cpu_efficiency"] = safe_ratio(cpu_time_sec, alloc_cpus * elapsed_sec)
        requested_mem_mb = row.get("requested_mem_mb")
        if requested_mem_mb:
            row["memory_efficiency"] = safe_ratio(row.get("max_rss_mb"), requested_mem_mb)
        requested_runtime_min = row.get("requested_runtime_min")
        if requested_runtime_min:
            row["runtime_efficiency"] = safe_ratio(elapsed_sec, requested_runtime_min * 60)
        state = str(row.get("state", "")).upper()
        row["oom"] = state.startswith("OUT_OF_MEMORY")
        row["timeout"] = state.startswith("TIMEOUT")
    return row

BENCHMARK_JOB_FIELDS = [
    "launch_index",
    "rule",
    "attempt",
    "logical_job_key",
    "internal_jobid",
    "external_jobid",
    "wildcards",
    "requested_cpus",
    "requested_mem_mb",
    "requested_runtime_min",
    "state",
    "exit_code",
    "alloc_cpus",
    "elapsed_sec",
    "cpu_time_sec",
    "max_rss_mb",
    "cpu_efficiency",
    "memory_efficiency",
    "runtime_efficiency",
    "oom",
    "timeout",
]

BENCHMARK_RULE_FIELDS = [
    "rule",
    "launches",
    "logical_jobs",
    "retries",
    "failed_launches",
    "oom_launches",
    "timeout_launches",
    "median_requested_cpus",
    "median_alloc_cpus",
    "median_requested_mem_mb",
    "median_max_rss_mb",
    "median_memory_efficiency",
    "median_requested_runtime_min",
    "median_elapsed_sec",
    "median_runtime_efficiency",
    "allocated_cpu_sec",
    "used_cpu_sec",
    "weighted_cpu_efficiency",
]

def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def write_benchmark_summary_file(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)

def summarize_benchmark_rows(job_rows):
    if not job_rows:
        return {
            "benchmarked_launches": 0,
            "logical_jobs": 0,
            "retries": 0,
            "failed_launches": 0,
            "oom_launches": 0,
            "timeout_launches": 0,
            "max_alloc_cpus": None,
            "peak_max_rss_mb": None,
            "total_elapsed_sec": 0,
            "allocated_cpu_sec": 0,
            "used_cpu_sec": 0,
            "weighted_cpu_efficiency": None,
            "jobs_missing_accounting": 0,
            "most_retried_rules": [],
        }

    attempt_rows = [row for row in job_rows if row.get("attempt", 1) > 1]
    peak_max_rss_mb = max((row["max_rss_mb"] for row in job_rows if row.get("max_rss_mb") is not None), default=None)
    max_alloc_cpus = max((row["alloc_cpus"] for row in job_rows if row.get("alloc_cpus") is not None), default=None)
    total_elapsed_sec = sum(row.get("elapsed_sec") or 0 for row in job_rows)
    allocated_cpu_sec = sum((row.get("alloc_cpus") or 0) * (row.get("elapsed_sec") or 0) for row in job_rows)
    used_cpu_sec = sum(row.get("cpu_time_sec") or 0 for row in job_rows)
    retry_counter = Counter(row["rule"] for row in attempt_rows)
    return {
        "benchmarked_launches": len(job_rows),
        "logical_jobs": len({row["logical_job_key"] for row in job_rows}),
        "retries": len(attempt_rows),
        "failed_launches": sum(1 for row in job_rows if str(row.get("state", "")).upper() not in {"", "COMPLETED"}),
        "oom_launches": sum(1 for row in job_rows if row.get("oom")),
        "timeout_launches": sum(1 for row in job_rows if row.get("timeout")),
        "max_alloc_cpus": max_alloc_cpus,
        "peak_max_rss_mb": peak_max_rss_mb,
        "total_elapsed_sec": total_elapsed_sec,
        "allocated_cpu_sec": allocated_cpu_sec,
        "used_cpu_sec": used_cpu_sec,
        "weighted_cpu_efficiency": safe_ratio(used_cpu_sec, allocated_cpu_sec),
        "jobs_missing_accounting": sum(1 for row in job_rows if row.get("state") is None),
        "most_retried_rules": retry_counter.most_common(5),
    }

def summarize_benchmark_rules(job_rows):
    grouped = defaultdict(list)
    for row in job_rows:
        grouped[row["rule"]].append(row)

    summaries = []
    for rule_name in sorted(grouped):
        rows = grouped[rule_name]
        summaries.append(
            {
                "rule": rule_name,
                "launches": len(rows),
                "logical_jobs": len({row["logical_job_key"] for row in rows}),
                "retries": sum(1 for row in rows if row.get("attempt", 1) > 1),
                "failed_launches": sum(1 for row in rows if str(row.get("state", "")).upper() not in {"", "COMPLETED"}),
                "oom_launches": sum(1 for row in rows if row.get("oom")),
                "timeout_launches": sum(1 for row in rows if row.get("timeout")),
                "median_requested_cpus": median_or_none(row.get("requested_cpus") for row in rows),
                "median_alloc_cpus": median_or_none(row.get("alloc_cpus") for row in rows),
                "median_requested_mem_mb": median_or_none(row.get("requested_mem_mb") for row in rows),
                "median_max_rss_mb": median_or_none(row.get("max_rss_mb") for row in rows),
                "median_memory_efficiency": median_or_none(row.get("memory_efficiency") for row in rows),
                "median_requested_runtime_min": median_or_none(row.get("requested_runtime_min") for row in rows),
                "median_elapsed_sec": median_or_none(row.get("elapsed_sec") for row in rows),
                "median_runtime_efficiency": median_or_none(row.get("runtime_efficiency") for row in rows),
                "allocated_cpu_sec": sum((row.get("alloc_cpus") or 0) * (row.get("elapsed_sec") or 0) for row in rows),
                "used_cpu_sec": sum(row.get("cpu_time_sec") or 0 for row in rows),
                "weighted_cpu_efficiency": safe_ratio(
                    sum(row.get("cpu_time_sec") or 0 for row in rows),
                    sum((row.get("alloc_cpus") or 0) * (row.get("elapsed_sec") or 0) for row in rows),
                ),
            }
        )
    return summaries

def write_benchmark_tables(paths, job_rows, rule_rows):
    write_tsv(paths["jobs"], BENCHMARK_JOB_FIELDS, job_rows)
    write_tsv(paths["rules"], BENCHMARK_RULE_FIELDS, rule_rows)

def generate_benchmark_reports(output_dir, run_id, job_rows):
    paths = build_benchmark_paths(output_dir, run_id)
    summary = summarize_benchmark_rows(job_rows)
    rule_rows = summarize_benchmark_rules(job_rows)

    write_benchmark_tables(paths, job_rows, rule_rows)
    write_benchmark_summary_file(
        paths["summary"],
        {
            "status": "generated",
            **summary,
            "allocated_cpu_hours": round(summary["allocated_cpu_sec"] / 3600, 4),
            "used_cpu_hours": round(summary["used_cpu_sec"] / 3600, 4),
            "rules": rule_rows,
        },
    )
    return {
        "paths": paths,
        "summary": summary,
        "rules": rule_rows,
    }

def generate_run_benchmark(output_dir, metadata_path=None, metadata=None, quiet=True):
    metadata = metadata or load_metadata_file(metadata_path)
    if not metadata:
        return None

    run_id = metadata.get("run_id")
    if not run_id:
        return None

    output_dir = Path(output_dir)
    paths = build_benchmark_paths(output_dir, run_id)
    base_summary = {
        "run_id": run_id,
        "command": metadata.get("command"),
        "profile": get_run_profile(metadata),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }

    if should_skip_benchmark(metadata):
        payload = {
            **base_summary,
            "status": "skipped",
            "message": "Resource benchmarking was skipped because --skip-benchmark was used for this run.",
        }
        write_benchmark_summary_file(paths["summary"], payload)
        update_launch_metadata(
            metadata_path,
            benchmark_status=payload["status"],
            benchmark_summary=str(paths["summary"]),
            benchmark_jobs=str(paths["jobs"]),
            benchmark_rules=str(paths["rules"]),
            benchmark_generated_at=payload["generated_at"],
        )
        return {
            "status": payload["status"],
            "paths": paths,
            "summary": None,
        }

    if not should_benchmark_run(metadata):
        payload = {
            **base_summary,
            "status": "unsupported_profile",
            "message": "Resource benchmarking currently supports only runs launched with the slurm profile.",
        }
        write_benchmark_summary_file(paths["summary"], payload)
        update_launch_metadata(
            metadata_path,
            benchmark_status=payload["status"],
            benchmark_summary=str(paths["summary"]),
            benchmark_jobs=str(paths["jobs"]),
            benchmark_rules=str(paths["rules"]),
            benchmark_generated_at=payload["generated_at"],
        )
        return {
            "status": payload["status"],
            "paths": paths,
            "summary": None,
        }

    configured_log = metadata.get("snakemake_log")
    if configured_log:
        snakemake_log_path = Path(configured_log)
    else:
        snakemake_log_path = build_snakemake_log_path(output_dir, run_id)
    if not snakemake_log_path.exists():
        payload = {
            **base_summary,
            "status": "log_missing",
            "message": "Snakemake log not found; no benchmark information could be generated.",
        }
        write_benchmark_summary_file(paths["summary"], payload)
        update_launch_metadata(
            metadata_path,
            benchmark_status=payload["status"],
            benchmark_summary=str(paths["summary"]),
            benchmark_jobs=str(paths["jobs"]),
            benchmark_rules=str(paths["rules"]),
            benchmark_generated_at=payload["generated_at"],
        )
        return {
            "status": payload["status"],
            "paths": paths,
            "summary": None,
        }

    launches = parse_snakemake_submitted_launches(snakemake_log_path)
    if not launches:
        write_benchmark_tables(paths, [], [])
        payload = {
            **base_summary,
            "status": "no_submitted_jobs",
            **summarize_benchmark_rows([]),
            "message": "No submitted non-local SLURM jobs were detected in the Snakemake log.",
            "rules": [],
        }
        write_benchmark_summary_file(paths["summary"], payload)
        result = {
            "status": payload["status"],
            "paths": paths,
            "summary": summarize_benchmark_rows([]),
        }
        update_launch_metadata(
            metadata_path,
            benchmark_status=result["status"],
            benchmark_summary=str(paths["summary"]),
            benchmark_jobs=str(paths["jobs"]),
            benchmark_rules=str(paths["rules"]),
            benchmark_generated_at=payload["generated_at"],
        )
        return result

    sacct_rows = query_sacct_for_jobs([launch["external_jobid"] for launch in launches])
    if sacct_rows is None:
        job_rows = [benchmark_job_row(launch, None) for launch in launches]
        rule_rows = summarize_benchmark_rules(job_rows)
        write_benchmark_tables(paths, job_rows, rule_rows)
        summary = summarize_benchmark_rows(job_rows)
        payload = {
            **base_summary,
            "status": "accounting_unavailable",
            **summary,
            "message": "SLURM sacct data could not be queried, so actual resource usage could not be summarized.",
            "rules": rule_rows,
        }
        write_benchmark_summary_file(paths["summary"], payload)
        result = {
            "status": payload["status"],
            "paths": paths,
            "summary": summary,
        }
        update_launch_metadata(
            metadata_path,
            benchmark_status=result["status"],
            benchmark_summary=str(paths["summary"]),
            benchmark_jobs=str(paths["jobs"]),
            benchmark_rules=str(paths["rules"]),
            benchmark_generated_at=payload["generated_at"],
        )
        if not quiet:
            print(f"{INFO}INFO:{RESET} SLURM accounting information is unavailable; benchmark files were not generated.")
        return result

    job_rows = [benchmark_job_row(launch, sacct_rows.get(launch["external_jobid"])) for launch in launches]
    result = generate_benchmark_reports(output_dir, run_id, job_rows)
    result["status"] = "generated"
    update_launch_metadata(
        metadata_path,
        benchmark_status=result["status"],
        benchmark_summary=str(paths["summary"]),
        benchmark_jobs=str(paths["jobs"]),
        benchmark_rules=str(paths["rules"]),
        benchmark_generated_at=base_summary["generated_at"],
    )
    return result
