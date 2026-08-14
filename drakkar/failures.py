from __future__ import annotations

import re
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path

from drakkar.benchmark import write_tsv
from drakkar.cli_context import ERROR, INFO, RESET
from drakkar.output import print, section
from drakkar.run_metadata import build_snakemake_log_path, load_metadata_file, update_launch_metadata

FAILURE_REPORT_FIELDS = [
    "run_id",
    "rule",
    "target",
    "attempts",
    "status",
    "category",
    "reason",
    "slurm_state",
    "internal_jobid",
    "external_jobid",
    "detail",
    "action",
    "job_log",
    "output",
    "last_failure_at",
]

# category -> (table label, what happened, what to do, resolution)
# resolution: "rerun" (relaunch as is), "scale" (relaunch with more resources), "fix" (needs changes)
FAILURE_CATEGORIES = {
    "timeout": (
        "timeout",
        "the job hit its SLURM wall-time limit",
        "Relaunch the same drakkar command with a larger --time-multiplier (e.g. --time-multiplier 2).",
        "scale",
    ),
    "out-of-memory": (
        "out-of-memory",
        "the job was killed after exceeding its memory allocation",
        "Relaunch the same drakkar command with a larger --memory-multiplier (e.g. --memory-multiplier 2).",
        "scale",
    ),
    "node-failure": (
        "node-failure",
        "the compute node failed or the job was preempted",
        "Transient cluster problem: relaunch the same drakkar command, no changes needed.",
        "rerun",
    ),
    "cancelled": (
        "cancelled",
        "the job was cancelled before it could finish",
        "Relaunch the same drakkar command once the cluster allows the job to run.",
        "rerun",
    ),
    "storage": (
        "storage",
        "the job ran out of disk space or exceeded a storage quota",
        "Free space or raise the quota on the output/scratch filesystem, then relaunch.",
        "fix",
    ),
    "missing-input": (
        "missing-input",
        "required input files were not found",
        "Check the sample information file, input directory, and database paths before relaunching.",
        "fix",
    ),
    "missing-output": (
        "missing-output",
        "the job finished but its output files were not visible",
        "Usually shared-filesystem latency: relaunch, and if it repeats raise --snakemake-latency-wait.",
        "rerun",
    ),
    "incomplete-files": (
        "incomplete",
        "output files from an interrupted run were left incomplete",
        "Relaunch with --snakemake-rerun-incomplete so Snakemake regenerates the incomplete files.",
        "rerun",
    ),
    "locked": (
        "locked",
        "the output directory is locked by a previous Snakemake run",
        "Run drakkar unlock -o <output_dir> and then relaunch.",
        "fix",
    ),
    "command-error": (
        "command-error",
        "the tool called by the rule exited with a non-zero status",
        "Inspect the job log of the failed rule and fix the cause before relaunching.",
        "fix",
    ),
    "unknown": (
        "unknown",
        "the failure could not be classified automatically",
        "Inspect the job log of the failed rule to find out what went wrong.",
        "fix",
    ),
}

RESOLUTION_ORDER = {"fix": 0, "scale": 1, "rerun": 2}

# Only the tail of a job log is inspected, so oversized logs stay cheap to read.
MAX_JOB_LOG_BYTES = 512 * 1024

ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*m")

ERROR_BLOCK_RE = re.compile(
    r"^Error (?:in|executing) rule\s+([A-Za-z_][A-Za-z0-9_]*)\b(.*)$"
)
BLOCK_START_RE = re.compile(r"^(localrule|localcheckpoint|rule|checkpoint)\s+(.+?):\s*$")
BLOCK_FIELD_RE = re.compile(r"^ {4}([A-Za-z][A-Za-z0-9_-]*):\s*(.*)$")
FINISHED_RE = re.compile(r"Finished\s+job(?:id:|\s+)\s*(\d+)")
SLURM_STATE_RE = re.compile(r"SLURM status is:\s*'?([A-Z_ ]+)'?")

ERROR_LINE_MARKERS = (
    "error",
    "exception",
    "traceback",
    "no such file",
    "not found",
    "permission denied",
    "cannot",
    "can't",
    "failed",
    "fatal",
    "aborted",
    "invalid",
    "out of memory",
    "segmentation fault",
    "core dumped",
    "no space left",
    "quota exceeded",
    "killed",
    "assertion",
    "unable to",
)

NOISE_LINE_MARKERS = (
    "check log file(s) for error details",
    "building dag of jobs",
    "select jobs to execute",
    "provided remote nodes",
    "using shell",
    "trying to restart",
    "will exit after finishing currently running jobs",
    "complete log:",
    "exiting because a job execution failed",
    "shutting down, this might take some time",
    "storing output in storage",
    "removing output files of failed job",
)

OOM_HINT_RE = re.compile(
    r"oom[-_ ]?kill|out of memory|exceeded job memory limit|memoryerror|std::bad_alloc|"
    r"cannot allocate memory|killed by signal 9",
    re.IGNORECASE,
)
TIMEOUT_HINT_RE = re.compile(r"due to time limit|time limit exceeded", re.IGNORECASE)
STORAGE_HINT_RE = re.compile(r"no space left on device|disk quota exceeded", re.IGNORECASE)

def build_failure_report_path(output_dir, run_id):
    return Path(output_dir) / "log" / f"drakkar_{run_id}.failures.tsv"

def clean_line(line):
    return ANSI_ESCAPE_RE.sub("", str(line or "")).replace("\r", "").rstrip("\n")

def truncate(text, width):
    text = " ".join(str(text or "").split())
    if len(text) <= width:
        return text
    if width <= 1:
        return text[:width]
    return text[: width - 1] + "…"

def parse_wildcard_values(text):
    values = []
    for item in str(text or "").split(","):
        if "=" not in item:
            continue
        values.append(item.split("=", 1)[1].strip())
    return [value for value in values if value]

def parse_resources(text):
    resources = {}
    for item in str(text or "").split(","):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        resources[key.strip()] = value.strip()
    return resources

def parse_log_paths(text):
    value = str(text or "").split(" (check log file")[0]
    return [item.strip() for item in value.split(",") if item.strip()]

def new_job_block(segment, rule):
    return {
        "segment": segment,
        "rule": rule,
        "jobid": None,
        "wildcards": "",
        "input": "",
        "output": "",
        "log": "",
        "resources": "",
    }

def is_noise_line(line):
    lowered = line.lower()
    return any(marker in lowered for marker in NOISE_LINE_MARKERS)

def looks_like_error_line(line):
    lowered = line.lower()
    return any(marker in lowered for marker in ERROR_LINE_MARKERS)

def parse_snakemake_failures(log_path):
    """Extract failed job events and workflow-level errors from a Snakemake log."""
    result = {"failures": [], "workflow_errors": [], "has_log": False}
    path = Path(log_path) if log_path else None
    if path is None or not path.exists():
        return result
    result["has_log"] = True

    segment = 0
    segment_has_activity = False
    current_block = None
    block_by_key = {}
    finished_at = {}
    failures = []
    workflow_errors = []
    current_error = None
    current_field = None
    recent_lines = deque(maxlen=40)
    last_timestamp = None
    pending_workflow_error = None
    line_count = 0

    def close_error(index):
        nonlocal current_error, current_field
        if current_error is not None:
            current_error["end_index"] = index
            failures.append(current_error)
        current_error = None
        current_field = None

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for index, raw_line in enumerate(handle):
            line_count = index + 1
            line = clean_line(raw_line)
            stripped = line.strip()

            if current_error is not None:
                if stripped and not line.startswith(" "):
                    close_error(index)
                else:
                    field_match = BLOCK_FIELD_RE.match(line)
                    if field_match:
                        current_field = field_match.group(1).lower()
                        current_error["fields"][current_field] = field_match.group(2).strip()
                    elif stripped and current_field:
                        existing = current_error["fields"].get(current_field, "")
                        if len(existing) < 4000:
                            current_error["fields"][current_field] = f"{existing}\n{stripped}".strip()
                    continue

            if pending_workflow_error is not None:
                keep_collecting = (
                    stripped
                    and pending_workflow_error["detail_lines"] < 6
                    and not is_noise_line(stripped)
                    and not ERROR_BLOCK_RE.match(stripped)
                    and not BLOCK_START_RE.match(stripped)
                )
                if keep_collecting:
                    pending_workflow_error["detail"] = f"{pending_workflow_error['detail']} {stripped}".strip()
                    pending_workflow_error["detail_lines"] += 1
                    recent_lines.append(stripped)
                    continue
                workflow_errors.append(pending_workflow_error)
                pending_workflow_error = None

            if not stripped:
                recent_lines.append(stripped)
                continue

            timestamp_match = re.fullmatch(r"\[(.+)\]", stripped)
            if timestamp_match:
                last_timestamp = timestamp_match.group(1)
                recent_lines.append(stripped)
                continue

            if stripped.startswith("Building DAG of jobs"):
                if segment_has_activity:
                    segment += 1
                    current_block = None
                    segment_has_activity = False
                recent_lines.append(stripped)
                continue

            error_match = ERROR_BLOCK_RE.match(stripped)
            if error_match:
                segment_has_activity = True
                trailer = error_match.group(2) or ""
                fields = {}
                jobid_match = re.search(r"jobid[:=]\s*(\d+)", trailer)
                if jobid_match:
                    fields["jobid"] = jobid_match.group(1)
                external_match = re.search(r"external(?:_jobid)?[:=]?\s*([0-9_]+)", trailer)
                if external_match:
                    fields["external_jobid"] = external_match.group(1)
                current_error = {
                    "rule": error_match.group(1),
                    "segment": segment,
                    "index": index,
                    "timestamp": last_timestamp,
                    "fields": fields,
                    "context": [item for item in recent_lines if item],
                }
                current_field = None
                recent_lines.append(stripped)
                continue

            block_match = BLOCK_START_RE.match(stripped)
            if block_match:
                segment_has_activity = True
                current_block = new_job_block(segment, block_match.group(2))
                recent_lines.append(stripped)
                continue

            if current_block is not None:
                field_match = BLOCK_FIELD_RE.match(line)
                if field_match:
                    field = field_match.group(1).lower()
                    value = field_match.group(2).strip()
                    if field == "jobid":
                        current_block["jobid"] = value
                        block_by_key[(segment, value)] = current_block
                    elif field in ("wildcards", "input", "output", "log", "resources"):
                        current_block[field] = value
                    recent_lines.append(stripped)
                    continue

            finished_match = FINISHED_RE.search(stripped)
            if finished_match:
                finished_at[(segment, finished_match.group(1))] = index
                recent_lines.append(stripped)
                continue

            workflow_error = classify_workflow_error_line(stripped)
            if workflow_error is not None:
                segment_has_activity = True
                pending_workflow_error = {
                    "segment": segment,
                    "index": index,
                    "category": workflow_error[0],
                    "title": stripped,
                    "detail": "",
                    "detail_lines": 0,
                    "rule": workflow_error[1],
                }

            recent_lines.append(stripped)

    close_error(line_count)
    if pending_workflow_error is not None:
        workflow_errors.append(pending_workflow_error)

    for failure in failures:
        fields = failure["fields"]
        jobid = fields.get("jobid")
        block = block_by_key.get((failure["segment"], jobid)) if jobid else None
        failure["internal_jobid"] = jobid or ""
        failure["external_jobid"] = str(fields.get("external_jobid", "")).strip()
        failure["message"] = fields.get("message", "")
        failure["output"] = fields.get("output", "") or (block or {}).get("output", "")
        failure["input"] = fields.get("input", "") or (block or {}).get("input", "")
        failure["wildcards"] = (block or {}).get("wildcards", "")
        failure["resources"] = parse_resources((block or {}).get("resources", ""))
        failure["job_logs"] = parse_log_paths(fields.get("log", "")) + parse_log_paths((block or {}).get("log", ""))
        failure["slurm_state"] = extract_slurm_state(failure["message"])
        failure["recovered"] = finished_at.get((failure["segment"], jobid), -1) > failure["index"] if jobid else False

    result["failures"] = failures
    result["workflow_errors"] = workflow_errors
    return result

def classify_workflow_error_line(line):
    """Return (category, rule) for workflow-level (non job) errors, or None."""
    rule_match = re.search(r"\brule\s+([A-Za-z_][A-Za-z0-9_]*)", line)
    rule = rule_match.group(1) if rule_match else ""
    if "MissingInputException" in line:
        return "missing-input", rule
    if "IncompleteFilesException" in line or "seem to be incomplete" in line:
        return "incomplete-files", rule
    if "LockException" in line or "Directory cannot be locked" in line:
        return "locked", rule
    if "MissingOutputException" in line:
        return "missing-output", rule
    if line.startswith("WorkflowError") or "InputFunctionException" in line or "ChildIOException" in line:
        return "unknown", rule
    return None

def extract_slurm_state(message):
    match = SLURM_STATE_RE.search(str(message or ""))
    if match:
        return match.group(1).strip()
    return ""

def read_job_log_details(job_logs, max_lines=400):
    """Return (detail, excerpt, hints) extracted from the SLURM/rule log of a failed job."""
    excerpt = []
    hints = set()
    for job_log in job_logs:
        path = Path(job_log)
        try:
            if not path.is_file():
                continue
            lines = deque(maxlen=max_lines)
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                size = path.stat().st_size
                if size > MAX_JOB_LOG_BYTES:
                    handle.seek(size - MAX_JOB_LOG_BYTES)
                    handle.readline()
                for raw_line in handle:
                    lines.append(clean_line(raw_line).strip())
        except OSError:
            continue

        text = "\n".join(lines)
        if OOM_HINT_RE.search(text):
            hints.add("out-of-memory")
        if TIMEOUT_HINT_RE.search(text):
            hints.add("timeout")
        if STORAGE_HINT_RE.search(text):
            hints.add("storage")

        candidates = [
            item for item in lines
            if item and looks_like_error_line(item) and not is_noise_line(item)
        ]
        if candidates:
            excerpt = candidates[-5:]
        elif not excerpt:
            excerpt = [item for item in lines if item][-3:]
    detail = excerpt[-1] if excerpt else ""
    return detail, excerpt, hints

def detail_from_context(context):
    for line in reversed(context or []):
        if not line or is_noise_line(line):
            continue
        if re.match(r"^(\[|rule |localrule |checkpoint |Job |Submitted |Select jobs|Execute |\d+ of \d+ steps)", line):
            continue
        if looks_like_error_line(line):
            return line
    return ""

def classify_failure(failure, hints=None):
    hints = hints or set()
    message = str(failure.get("message", ""))
    state = str(failure.get("slurm_state", "")).upper()
    shell = str(failure.get("fields", {}).get("shell", ""))
    combined = f"{message}\n{shell}"

    if "out-of-memory" in hints or state.startswith("OUT_OF_MEMORY") or "OUT_OF_ME" in state:
        return "out-of-memory"
    if state.startswith("TIMEOUT") or "timeout" in hints or "DUE TO TIME LIMIT" in message.upper():
        return "timeout"
    if "storage" in hints:
        return "storage"
    if state.startswith("NODE_FAIL") or state.startswith("PREEMPTED") or state.startswith("BOOT_FAIL"):
        return "node-failure"
    if state.startswith("CANCELLED"):
        return "cancelled"
    if "MissingInputException" in combined:
        return "missing-input"
    if "MissingOutputException" in combined or "Missing files after" in combined:
        return "missing-output"
    if "IncompleteFilesException" in combined:
        return "incomplete-files"
    if "LockException" in combined:
        return "locked"
    if "non-zero exit code" in combined or state.startswith("FAILED") or failure.get("external_jobid"):
        return "command-error"
    if message or shell:
        return "command-error"
    return "unknown"

def resource_detail(failure, category):
    """Describe resource-driven failures with the limit the job actually requested."""
    resources = failure.get("resources") or {}
    state = str(failure.get("slurm_state", "")).upper() or "unknown"
    if category == "timeout":
        runtime = resources.get("runtime")
        if runtime:
            return f"hit the SLURM time limit ({runtime} min requested)"
        return "hit the SLURM time limit"
    if category == "out-of-memory":
        memory = resources.get("mem_mb") or resources.get("mem_mib")
        if memory:
            return f"exceeded its memory allocation ({memory} MB requested)"
        return "exceeded its memory allocation"
    if category in ("node-failure", "cancelled"):
        return f"SLURM reported state {state}"
    return ""

def failure_target(failure):
    values = parse_wildcard_values(failure.get("wildcards", ""))
    if values:
        return ",".join(values)
    output = str(failure.get("output", "")).split(",")[0].strip()
    if output:
        return Path(output).name
    jobid = failure.get("internal_jobid")
    return f"job {jobid}" if jobid else "-"

def build_failure_groups(failures, read_job_logs=True):
    groups = {}
    for failure in failures:
        target = failure_target(failure)
        key = (failure.get("rule"), target)
        group = groups.get(key)
        if group is None:
            group = {
                "rule": failure.get("rule"),
                "target": target,
                "attempts": 0,
                "events": [],
            }
            groups[key] = group
        group["attempts"] += 1
        group["events"].append(failure)

    rows = []
    for group in groups.values():
        last = group["events"][-1]
        detail, excerpt, hints = ("", [], set())
        if read_job_logs:
            detail, excerpt, hints = read_job_log_details(last.get("job_logs") or [])
        category = classify_failure(last, hints)
        resource_text = resource_detail(last, category)
        if resource_text:
            detail = resource_text
        if not detail:
            detail = detail_from_context(last.get("context"))
        if not detail:
            detail = str(last.get("message", "")).split("For further error details")[0].strip()
        if not excerpt and detail:
            excerpt = [detail]
        label, reason, action, resolution = FAILURE_CATEGORIES.get(category, FAILURE_CATEGORIES["unknown"])
        rows.append(
            {
                "rule": group["rule"],
                "target": group["target"],
                "attempts": group["attempts"],
                "status": "recovered" if last.get("recovered") else "failed",
                "category": category,
                "label": label,
                "reason": reason,
                "action": action,
                "resolution": resolution,
                "detail": detail,
                "excerpt": excerpt,
                "slurm_state": last.get("slurm_state", ""),
                "internal_jobid": last.get("internal_jobid", ""),
                "external_jobid": last.get("external_jobid", ""),
                "job_logs": last.get("job_logs") or [],
                "output": last.get("output", ""),
                "last_failure_at": last.get("timestamp") or "",
            }
        )

    rows.sort(
        key=lambda row: (
            row["status"] != "failed",
            RESOLUTION_ORDER.get(row["resolution"], 9),
            row["rule"] or "",
            row["target"] or "",
        )
    )
    return rows

def summarize_failure_rows(rows, workflow_errors=None):
    failed = [row for row in rows if row["status"] == "failed"]
    recovered = [row for row in rows if row["status"] == "recovered"]
    categories = Counter(row["category"] for row in failed)
    resolutions = {row["resolution"] for row in failed}
    for error in workflow_errors or []:
        categories[error.get("category", "unknown")] += 0
    if not failed:
        verdict = "none"
    elif "fix" in resolutions:
        verdict = "needs-changes"
    elif "scale" in resolutions:
        verdict = "rerun-with-more-resources"
    else:
        verdict = "rerun"
    return {
        "failed_jobs": len(failed),
        "failed_rules": len({row["rule"] for row in failed}),
        "recovered_jobs": len(recovered),
        "total_attempts": sum(row["attempts"] for row in rows),
        "categories": dict(categories),
        "verdict": verdict,
    }

def build_failure_report(log_path, run_id=None, read_job_logs=True):
    parsed = parse_snakemake_failures(log_path)
    rows = build_failure_groups(parsed["failures"], read_job_logs=read_job_logs)
    workflow_errors = parsed["workflow_errors"]
    return {
        "run_id": run_id,
        "has_log": parsed["has_log"],
        "log_path": str(log_path) if log_path else None,
        "rows": rows,
        "workflow_errors": workflow_errors,
        "summary": summarize_failure_rows(rows, workflow_errors),
        "report_path": None,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }

def failure_report_rows_for_tsv(report):
    run_id = report.get("run_id") or ""
    rows = []
    for row in report["rows"]:
        rows.append(
            {
                "run_id": run_id,
                "rule": row["rule"],
                "target": row["target"],
                "attempts": row["attempts"],
                "status": row["status"],
                "category": row["category"],
                "reason": row["reason"],
                "slurm_state": row["slurm_state"],
                "internal_jobid": row["internal_jobid"],
                "external_jobid": row["external_jobid"],
                "detail": truncate(row["detail"], 300),
                "action": row["action"],
                "job_log": row["job_logs"][0] if row["job_logs"] else "",
                "output": row["output"],
                "last_failure_at": row["last_failure_at"],
            }
        )
    for error in report.get("workflow_errors") or []:
        label, reason, action, _resolution = FAILURE_CATEGORIES.get(
            error.get("category", "unknown"), FAILURE_CATEGORIES["unknown"]
        )
        rows.append(
            {
                "run_id": run_id,
                "rule": error.get("rule") or "-",
                "target": "workflow",
                "attempts": 1,
                "status": "failed",
                "category": error.get("category", "unknown"),
                "reason": reason,
                "slurm_state": "",
                "internal_jobid": "",
                "external_jobid": "",
                "detail": truncate(f"{error.get('title', '')} {error.get('detail', '')}", 300),
                "action": action,
                "job_log": "",
                "output": "",
                "last_failure_at": "",
            }
        )
    return rows

def print_failure_report(report, detailed=False, limit=15):
    section("FAILURE REPORT")
    if not report.get("has_log"):
        print("No Snakemake log was found, so no failure report could be built.")
        return
    rows = report["rows"]
    workflow_errors = report.get("workflow_errors") or []
    summary = report["summary"]
    if not rows and not workflow_errors:
        print("No failed jobs were detected in the Snakemake log.")
        return

    failed_rows = [row for row in rows if row["status"] == "failed"]
    print(
        f"Failed jobs: {summary['failed_jobs']} across {summary['failed_rules']} rule(s) "
        f"({summary['total_attempts']} failed attempts in total)"
    )
    if summary["recovered_jobs"]:
        print(f"Failed once but completed after a retry: {summary['recovered_jobs']}")

    if failed_rows:
        print("")
        print(f"{'Rule':<26} {'Target':<18} {'Att':>3}  {'Reason':<14}  Detail")
        shown = failed_rows if detailed or limit is None else failed_rows[:limit]
        for row in shown:
            print(
                f"{truncate(row['rule'], 26):<26} {truncate(row['target'], 18):<18} "
                f"{row['attempts']:>3}  {truncate(row['label'], 14):<14}  {truncate(row['detail'], 60)}"
            )
        hidden = len(failed_rows) - len(shown)
        if hidden > 0:
            print(f"... and {hidden} more failed job(s). Use --failures to list them all.")

    if workflow_errors:
        print("")
        print("Workflow-level errors:")
        for error in workflow_errors:
            label = FAILURE_CATEGORIES.get(error.get("category", "unknown"), FAILURE_CATEGORIES["unknown"])[0]
            print(f"  [{label}] {truncate(error.get('title', ''), 110)}")
            if error.get("detail"):
                print(f"      {truncate(error['detail'], 110)}")

    print_failure_actions(report)

    if detailed:
        print_failure_details(report)

    report_path = report.get("report_path")
    if report_path:
        print("")
        print(f"Full failure table: {report_path}")

def print_failure_actions(report):
    rows = [row for row in report["rows"] if row["status"] == "failed"]
    workflow_errors = report.get("workflow_errors") or []
    if not rows and not workflow_errors:
        return

    grouped = {}
    for row in rows:
        grouped.setdefault(row["category"], []).append(row)
    for error in workflow_errors:
        grouped.setdefault(error.get("category", "unknown"), [])

    section("WHAT TO DO NEXT")
    for category in sorted(
        grouped,
        key=lambda name: RESOLUTION_ORDER.get(FAILURE_CATEGORIES.get(name, FAILURE_CATEGORIES["unknown"])[3], 9),
    ):
        category_rows = grouped[category]
        label, reason, action, _resolution = FAILURE_CATEGORIES.get(category, FAILURE_CATEGORIES["unknown"])
        rule_names = sorted({row["rule"] for row in category_rows if row.get("rule")})
        rules_text = ", ".join(rule_names[:4])
        if len(rule_names) > 4:
            rules_text += f", +{len(rule_names) - 4} more"
        if category_rows:
            print(f"[{label}] {len(category_rows)} job(s) in {rules_text or 'unknown rule'}: {reason}.")
        else:
            print(f"[{label}] {reason}.")
        print(f"  -> {action}")
        job_log = next((row["job_logs"][0] for row in category_rows if row.get("job_logs")), None)
        if job_log and FAILURE_CATEGORIES.get(category, FAILURE_CATEGORIES["unknown"])[3] == "fix":
            print(f"  -> Example job log: {job_log}")

    verdict = report["summary"]["verdict"]
    print("")
    if verdict == "rerun":
        print("Verdict: relaunching the same drakkar command should be enough.")
    elif verdict == "rerun-with-more-resources":
        print("Verdict: relaunch the same drakkar command with the resource multipliers suggested above.")
    elif verdict == "needs-changes":
        print("Verdict: at least one failure needs to be inspected and fixed before relaunching.")
    else:
        print("Verdict: no unresolved failures were detected.")
    print("Completed outputs are kept, so a relaunch resumes from where the workflow stopped.")

def print_failure_details(report):
    rows = [row for row in report["rows"] if row["status"] == "failed"]
    if not rows:
        return
    section("FAILED JOB DETAILS")
    for row in rows:
        print(f"{row['rule']} [{row['target']}] - {row['label']} ({row['attempts']} attempt(s))")
        if row["slurm_state"]:
            print(f"  SLURM state: {row['slurm_state']} (job {row['external_jobid'] or 'unknown'})")
        if row["output"]:
            print(f"  Expected output: {truncate(row['output'], 160)}")
        for job_log in row["job_logs"][:2]:
            print(f"  Job log: {job_log}")
        for line in row["excerpt"]:
            print(f"    | {truncate(line, 160)}")
        print("")

def generate_failure_report(output_dir, metadata_path=None, metadata=None, write=True, read_job_logs=True):
    metadata = metadata or load_metadata_file(metadata_path)
    if not metadata:
        return None
    run_id = metadata.get("run_id")
    if not run_id:
        return None

    output_dir = Path(output_dir)
    configured_log = metadata.get("snakemake_log")
    log_path = Path(configured_log) if configured_log else build_snakemake_log_path(output_dir, run_id)
    report = build_failure_report(log_path, run_id=run_id, read_job_logs=read_job_logs)

    if write and (report["rows"] or report["workflow_errors"]):
        report_path = build_failure_report_path(output_dir, run_id)
        try:
            write_tsv(report_path, FAILURE_REPORT_FIELDS, failure_report_rows_for_tsv(report))
            report["report_path"] = str(report_path)
        except OSError:
            report["report_path"] = None
        if metadata_path:
            update_launch_metadata(
                metadata_path,
                failure_report=report["report_path"],
                failed_jobs=report["summary"]["failed_jobs"],
                failure_verdict=report["summary"]["verdict"],
            )
    return report

def report_run_failures(output_dir, metadata_path=None, metadata=None, detailed=False):
    """Build, store, and print the failure report for a finished run."""
    try:
        report = generate_failure_report(output_dir, metadata_path=metadata_path, metadata=metadata)
    except Exception as exc:  # pragma: no cover - report generation must never mask the run failure
        print(f"{INFO}INFO:{RESET} Failure report could not be generated ({exc.__class__.__name__}: {exc}).")
        return None
    if report is None:
        return None
    if not report["rows"] and not report["workflow_errors"]:
        return report
    print_failure_report(report, detailed=detailed)
    return report
