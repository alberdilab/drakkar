import re
import shlex
from collections import Counter, deque
from pathlib import Path

import yaml

from drakkar.benchmark import format_hours, format_megabytes, format_percent, generate_run_benchmark
from drakkar.cli_context import ERROR, INFO, RESET, WORKFLOW_RUN_COMMANDS
from drakkar.output import print, section
from drakkar.run_metadata import build_snakemake_log_path, load_metadata_file
from drakkar.system_checks import is_snakemake_locked

def workflow_run_sort_key(metadata_path):
    path = Path(metadata_path)
    match = re.search(r"drakkar_(\d{8}-\d{6})\.ya?ml$", path.name)
    if match:
        return match.group(1)
    return path.name

def is_launch_metadata_path(metadata_path):
    return re.fullmatch(r"drakkar_\d{8}-\d{6}\.ya?ml", Path(metadata_path).name) is not None

def run_id_from_metadata_name(value):
    name = Path(str(value)).name
    match = re.fullmatch(r"drakkar_(\d{8}-\d{6})(?:_resources)?\.ya?ml", name)
    if match:
        return match.group(1)
    return str(value).strip().removeprefix("drakkar_").removesuffix(".yaml").removesuffix(".yml")

def discover_run_metadata(output_dir):
    output_path = Path(output_dir)
    runs = []
    for metadata_path in sorted(output_path.glob("drakkar_*.yaml"), key=workflow_run_sort_key, reverse=True):
        if not is_launch_metadata_path(metadata_path):
            continue
        try:
            with open(metadata_path, "r", encoding="utf-8") as handle:
                metadata = yaml.safe_load(handle) or {}
        except OSError:
            continue
        command = metadata.get("command")
        if command not in WORKFLOW_RUN_COMMANDS:
            continue
        runs.append((metadata_path, metadata))
    return runs

def resolve_run_metadata(output_dir, run_id=None):
    runs = discover_run_metadata(output_dir)
    if not runs:
        return None, None
    if not run_id:
        return runs[0]

    selector = str(run_id).strip()
    normalized_selector = run_id_from_metadata_name(selector)
    for metadata_path, metadata in runs:
        run_value = str(metadata.get("run_id", "")).strip()
        if selector == metadata_path.name or normalized_selector == run_value:
            return metadata_path, metadata
    return None, None

def discover_snakemake_fallback_logs(output_dir):
    log_dir = Path(output_dir) / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    return sorted((path for path in log_dir.glob("*") if path.is_file()), key=lambda path: path.stat().st_mtime, reverse=True)

def tail_file(path, line_count):
    lines = deque(maxlen=line_count)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lines.append(line.rstrip("\n"))
    return list(lines)

def extract_failure_excerpt(path, line_count=40):
    markers = (
        "RuleException",
        "MissingInputException",
        "WorkflowError",
        "CalledProcessError",
        "LockException",
        "Error in rule",
        "Traceback (most recent call last):",
    )
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle]
    last_index = None
    for index, line in enumerate(lines):
        if any(marker in line for marker in markers):
            last_index = index
    if last_index is None:
        return []
    start = max(0, last_index - 5)
    end = min(len(lines), last_index + line_count)
    return lines[start:end]

def classify_error_line(line):
    known_error_types = (
        "RuleException",
        "MissingInputException",
        "WorkflowError",
        "CalledProcessError",
        "LockException",
        "InputFunctionException",
        "ChildIOException",
    )
    if line.startswith("Error in rule "):
        return "RuleError"
    for error_type in known_error_types:
        if re.search(rf"\b{re.escape(error_type)}\b", line):
            return error_type
    match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*(?:Exception|Error))(?::|\b)", line)
    if match:
        return match.group(1)
    return None

def summarize_snakemake_log(path, metadata=None):
    summary = {
        "planned_jobs": None,
        "completed_steps": None,
        "total_steps": None,
        "percent_complete": None,
        "unique_rules": 0,
        "rule_executions": 0,
        "failed_rules": 0,
        "error_types": Counter(),
        "most_active_rules": [],
    }
    if path is None or not Path(path).exists():
        return summary

    rule_counter = Counter()
    finished_job_ids = set()
    started_job_ids = set()
    planned_jobs = None
    best_progress = None
    failed_rules = 0
    error_types = Counter()

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            match = re.match(r"^(?:local)?rule\s+(.+?):\s*$", stripped)
            if match:
                rule_counter[match.group(1)] += 1

            match = re.search(r"\bjobid:\s*(\d+)\b", stripped)
            if match:
                started_job_ids.add(match.group(1))

            match = re.search(r"Finished\s+job(?:id:|\s+)\s*(\d+)", stripped)
            if match:
                finished_job_ids.add(match.group(1))

            match = re.match(r"(\d+)\s+of\s+(\d+)\s+steps\s+\((\d+)%\)\s+done", stripped)
            if match:
                progress = (int(match.group(1)), int(match.group(2)), int(match.group(3)))
                if best_progress is None or progress[2] > best_progress[2] or (
                    progress[2] == best_progress[2] and progress[0] > best_progress[0]
                ):
                    best_progress = progress

            match = re.match(r"total\s+(\d+)\s*$", stripped)
            if match:
                count = int(match.group(1))
                planned_jobs = count if planned_jobs is None else max(planned_jobs, count)

            if stripped.startswith("Error in rule "):
                failed_rules += 1

            error_type = classify_error_line(stripped)
            if error_type:
                error_types[error_type] += 1

    completed_steps = best_progress[0] if best_progress else None
    total_steps = best_progress[1] if best_progress else None
    percent_complete = best_progress[2] if best_progress else None

    if completed_steps is None and finished_job_ids:
        completed_steps = len(finished_job_ids)
    if total_steps is None and planned_jobs is not None:
        total_steps = planned_jobs
    if percent_complete is None and completed_steps is not None and total_steps:
        percent_complete = int(round((completed_steps / total_steps) * 100))

    status = str((metadata or {}).get("status", "")).strip().lower()
    if status == "success":
        if total_steps is None:
            total_steps = planned_jobs or len(finished_job_ids) or len(started_job_ids) or None
        if total_steps is not None:
            completed_steps = total_steps
            percent_complete = 100

    summary.update(
        {
            "planned_jobs": planned_jobs,
            "completed_steps": completed_steps,
            "total_steps": total_steps,
            "percent_complete": percent_complete,
            "unique_rules": len(rule_counter),
            "rule_executions": sum(rule_counter.values()),
            "failed_rules": failed_rules,
            "error_types": error_types,
            "most_active_rules": rule_counter.most_common(5),
        }
    )
    return summary

def print_snakemake_summary(summary):
    section("EXECUTION SUMMARY")
    planned_jobs = summary.get("planned_jobs")
    completed_steps = summary.get("completed_steps")
    total_steps = summary.get("total_steps")
    percent_complete = summary.get("percent_complete")
    unique_rules = summary.get("unique_rules", 0)
    rule_executions = summary.get("rule_executions", 0)
    failed_rules = summary.get("failed_rules", 0)
    error_types = summary.get("error_types") or Counter()
    most_active_rules = summary.get("most_active_rules") or []

    print(f"Planned jobs: {planned_jobs if planned_jobs is not None else 'unknown'}")
    if percent_complete is not None and completed_steps is not None and total_steps is not None:
        print(f"Workflow progress: {percent_complete}% ({completed_steps}/{total_steps} steps)")
    elif completed_steps is not None:
        print(f"Completed steps observed: {completed_steps}")
    else:
        print("Workflow progress: unknown")
    print(f"Rules observed: {unique_rules} unique, {rule_executions} executions")
    print(f"Failed rules detected: {failed_rules}")
    if error_types:
        formatted_types = ", ".join(f"{name} ({count})" for name, count in error_types.most_common())
        print(f"Error types: {formatted_types}")
    else:
        print("Error types: none detected")
    if most_active_rules:
        formatted_rules = ", ".join(f"{name} ({count})" for name, count in most_active_rules)
        print(f"Most active rules: {formatted_rules}")

def print_benchmark_summary(benchmark_result):
    if not benchmark_result:
        return

    status = benchmark_result.get("status")
    if status == "skipped":
        section("RESOURCE BENCHMARK")
        print("Resource benchmarking was skipped for this run because --skip-benchmark was used.")
        return
    if status == "accounting_unavailable":
        section("RESOURCE BENCHMARK")
        print("SLURM accounting data could not be queried, so actual usage metrics are unavailable for this run.")
        return
    if status == "unsupported_profile":
        section("RESOURCE BENCHMARK")
        print("Resource benchmarking is currently available only for runs launched with the slurm profile.")
        return
    if status == "log_missing":
        section("RESOURCE BENCHMARK")
        print("Snakemake log not found, so no resource benchmark summary is available for this run.")
        return

    summary = benchmark_result.get("summary") or {}
    if status == "no_submitted_jobs":
        section("RESOURCE BENCHMARK")
        print("No submitted non-local SLURM jobs were detected in the Snakemake log.")
        return

    section("RESOURCE BENCHMARK")
    print(f"Benchmarked launches: {summary.get('benchmarked_launches', 0)}")
    print(f"Logical jobs: {summary.get('logical_jobs', 0)}")
    print(f"Relaunches detected: {summary.get('retries', 0)}")
    print(f"Failed launches: {summary.get('failed_launches', 0)}")
    print(
        f"Timeouts / OOM: {summary.get('timeout_launches', 0)} / {summary.get('oom_launches', 0)}"
    )
    print(f"Jobs missing sacct data: {summary.get('jobs_missing_accounting', 0)}")
    print(f"Peak allocated CPUs (single launch): {summary.get('max_alloc_cpus', 'unknown')}")
    print(f"Peak memory observed: {format_megabytes(summary.get('peak_max_rss_mb'))}")
    print(f"Total elapsed wall time: {format_hours((summary.get('total_elapsed_sec') or 0) / 3600)} hours")
    print(
        "Allocated CPU time: "
        f"{format_hours((summary.get('allocated_cpu_sec') or 0) / 3600)} CPU-hours"
    )
    print(
        "Used CPU time: "
        f"{format_hours((summary.get('used_cpu_sec') or 0) / 3600)} CPU-hours"
    )
    print(f"Weighted CPU efficiency: {format_percent(summary.get('weighted_cpu_efficiency'))}")
    most_retried = summary.get("most_retried_rules") or []
    if most_retried:
        formatted = ", ".join(f"{name} ({count})" for name, count in most_retried)
        print(f"Most retried rules: {formatted}")

def print_logging_usage_guide(output_path, selected_run_id=None):
    section("HOW TO INSPECT MORE")
    output_arg = shlex.quote(str(output_path))
    summary_cmd = f"drakkar logging -o {output_arg} --summary"
    paths_cmd = f"drakkar logging -o {output_arg} --paths"
    excerpt_cmd = f"drakkar logging -o {output_arg} --excerpt"
    full_cmd = f"drakkar logging -o {output_arg} --full"
    list_cmd = f"drakkar logging -o {output_arg} --list"
    print("Use these commands to inspect more detail:")
    print(f"  Summary only: {summary_cmd}")
    if selected_run_id:
        run_cmd = f"drakkar logging -o {output_arg} --run {selected_run_id} --summary"
        print(f"  Specific run summary: {run_cmd}")
    print(f"  Relevant paths: {paths_cmd}")
    print(f"  Failure excerpt or tail: {excerpt_cmd}")
    print(f"  Full Snakemake log: {full_cmd}")
    print(f"  Available runs: {list_cmd}")
    print("  Benchmark tables: inspect benchmark/drakkar_<run_id>.jobs.tsv and .rules.tsv when present")

def run_logging(output_dir, run_id=None, tail=50, full=False, paths=False, list_runs=False, summary=False, excerpt=False):
    output_path = Path(output_dir).resolve()
    if not output_path.exists():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path is not a directory: {output_path}")
        return 1

    section("DRAKKAR LOGGING")
    print(f"Output directory: {output_path}")
    print(f"Locked: {'yes' if is_snakemake_locked(str(output_path)) else 'no'}")

    runs = discover_run_metadata(output_path)
    if list_runs:
        if not runs:
            print(f"{INFO}INFO:{RESET} No workflow run metadata found in {output_path}.")
        else:
            section("AVAILABLE RUNS")
            for metadata_path, metadata in runs:
                run_value = metadata.get("run_id", metadata_path.stem.removeprefix("drakkar_"))
                command = metadata.get("command", "unknown")
                status = metadata.get("status", "unknown")
                print(f"{run_value}: command={command}, status={status}")
        return 0

    metadata_path, metadata = resolve_run_metadata(output_path, run_id)
    if run_id and metadata is None:
        print(f"{ERROR}ERROR:{RESET} Run not found in {output_path}: {run_id}")
        return 1
    snakemake_log_path = None
    fallback_logs = discover_snakemake_fallback_logs(output_path)
    if metadata is not None:
        configured_log = metadata.get("snakemake_log")
        if configured_log:
            snakemake_log_path = Path(configured_log)
        elif metadata.get("run_id"):
            candidate = build_snakemake_log_path(output_path, metadata["run_id"])
            if candidate.exists():
                snakemake_log_path = candidate
    if snakemake_log_path is None and fallback_logs:
        snakemake_log_path = fallback_logs[0]

    if metadata is not None:
        benchmark_result = generate_run_benchmark(output_path, metadata_path=metadata_path, metadata=metadata, quiet=True)
        metadata = load_metadata_file(metadata_path) or metadata
        section("RUN SUMMARY")
        selected_run_id = metadata.get('run_id', metadata_path.stem.removeprefix('drakkar_'))
        print(f"Run ID: {selected_run_id}")
        print(f"Command: {metadata.get('command', 'unknown')}")
        modules = metadata.get("modules") or []
        print(f"Modules: {', '.join(modules) if modules else 'unknown'}")
        print(f"Status: {metadata.get('status', 'unknown')}")
        print(f"Started: {metadata.get('started_at', metadata.get('timestamp', 'unknown'))}")
        if metadata.get("finished_at"):
            print(f"Finished: {metadata['finished_at']}")
        if "exit_code" in metadata:
            print(f"Exit code: {metadata['exit_code']}")
        if metadata.get("current_workflow"):
            print(f"Current workflow: {metadata['current_workflow']}")
        print(f"Metadata file: {metadata_path}")
    else:
        selected_run_id = None
        benchmark_result = None
        print(f"{INFO}INFO:{RESET} No workflow metadata found. Falling back to Snakemake logs only.")

    log_summary = summarize_snakemake_log(snakemake_log_path, metadata=metadata)
    if snakemake_log_path is not None and snakemake_log_path.exists():
        print_snakemake_summary(log_summary)
    if benchmark_result is not None:
        print_benchmark_summary(benchmark_result)

    if paths:
        section("LOG PATHS")
        if snakemake_log_path is not None:
            print(f"Main Snakemake log: {snakemake_log_path}")
        else:
            print("Main Snakemake log: not found")
        for fallback_log in fallback_logs[:5]:
            print(f"Snakemake fallback log: {fallback_log}")
        extra_logs = sorted(
            path for path in (output_path / "log").rglob("*")
            if path.is_file() and path != snakemake_log_path
        ) if (output_path / "log").exists() else []
        for extra_log in extra_logs[:20]:
            print(f"Additional log: {extra_log}")
        if metadata is not None:
            for key in ("benchmark_jobs", "benchmark_rules", "benchmark_summary"):
                benchmark_path = metadata.get(key)
                if benchmark_path:
                    print(f"{key.replace('_', ' ').title()}: {benchmark_path}")

    if snakemake_log_path is None or not snakemake_log_path.exists():
        print(f"{INFO}INFO:{RESET} No Snakemake log file found in {output_path}.")
        return 0

    if not full and not excerpt:
        print_logging_usage_guide(output_path, selected_run_id=selected_run_id)
        if summary:
            return 0
        if not paths:
            return 0

    if summary and not full and not excerpt:
        return 0

    section("SNAKEMAKE LOG")
    print(f"Log file: {snakemake_log_path}")
    if full:
        with open(snakemake_log_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                print(line, end="")
        return 0

    excerpt = extract_failure_excerpt(snakemake_log_path)
    if excerpt:
        print("Most recent failure excerpt:")
        for line in excerpt:
            print(line)
        return 0

    print(f"Last {tail} lines:")
    for line in tail_file(snakemake_log_path, tail):
        print(line)
    return 0
