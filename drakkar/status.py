from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from drakkar.cli_context import ERROR, INFO, PACKAGE_DIR, RESET
from drakkar.output import print, section
from drakkar.run_logs import (
    discover_snakemake_fallback_logs,
    resolve_run_metadata,
)
from drakkar.run_metadata import build_snakemake_log_path, load_metadata_file


HELPER_RULES = {
    "all",
    "create_report",
    "concatenate_or_link_reads",
    "concatenate_or_link_preprocessed",
    "samtools_stats",
    "assembly_flagstat",
    "maxbin2_table",
    "semibin2_table",
    "comebin_table",
    "rename_bins",
    "move_metadata",
    "rename_derep_headers",
    "profiling_total_reads",
    "split_coverm",
    "singlem_merge",
    "checkm2_metadata",
    "gtdbtk_input",
    "select_kegg",
    "dbcan2",
    "dbcan3",
    "dbcan4",
    "dbcan_summary",
    "antismash_regions",
    "genomad_regions",
    "merge_gene_annotations",
    "merge_cluster_annotations",
    "final_gene_annotation_table",
    "final_cluster_annotation_table",
    "merge_metagenome_gff",
    "write_database_versions",
}


WORKFLOW_RULE_FILES = {
    "preprocessing": [
        "rules/preparing.smk",
        "rules/preprocessing_ref.smk",
        "rules/preprocessing.smk",
        "rules/preprocessing_nonpareil.smk",
    ],
    "cataloging": ["rules/cataloging.smk"],
    "profiling": [
        "rules/profiling_genomes.smk",
        "rules/profiling_pangenomes.smk",
    ],
    "dereplicating": ["rules/dereplicating.smk"],
    "annotating": [
        "rules/annotating_taxonomy.smk",
        "rules/annotating_function.smk",
        "rules/annotating_network.smk",
    ],
    "inspecting": ["rules/inspecting_microdiversity.smk"],
    "expressing": ["rules/expressing.smk"],
    "database": ["rules/databases.smk"],
    "environments": ["rules/environments.smk"],
}


def _looks_like_metadata_file(value):
    text = str(value or "").strip()
    if not text:
        return False
    path = Path(text)
    return path.suffix.lower() in {".yaml", ".yml"} or re.match(
        r"^drakkar_\d{8}-\d{6}\.ya?ml$", path.name
    )


def resolve_status_target(target=None, output_dir=None, run_id=None):
    output_path = Path(output_dir).expanduser() if output_dir else Path.cwd()
    run_selector = run_id

    if target:
        target_path = Path(target).expanduser()
        if _looks_like_metadata_file(target):
            metadata_path = target_path
            if not metadata_path.is_absolute():
                cwd_candidate = Path.cwd() / metadata_path
                output_candidate = output_path / metadata_path
                metadata_path = cwd_candidate if cwd_candidate.exists() else output_candidate
            output_path = metadata_path.parent
            run_selector = metadata_path.name
        else:
            output_path = target_path
    elif run_id and _looks_like_metadata_file(run_id) and (
        Path(str(run_id)).is_absolute() or "/" in str(run_id)
    ):
        metadata_path = Path(str(run_id)).expanduser()
        if not metadata_path.is_absolute():
            cwd_candidate = Path.cwd() / metadata_path
            output_candidate = output_path / metadata_path
            metadata_path = cwd_candidate if cwd_candidate.exists() else output_candidate
        output_path = metadata_path.parent
        run_selector = metadata_path.name

    return output_path.resolve(), run_selector


def resolve_status_metadata(output_path, run_selector=None):
    if run_selector and (Path(str(run_selector)).is_absolute() or "/" in str(run_selector)):
        metadata_path = Path(run_selector).expanduser()
        if not metadata_path.is_absolute():
            cwd_candidate = Path.cwd() / metadata_path
            output_candidate = Path(output_path) / metadata_path
            metadata_path = cwd_candidate if cwd_candidate.exists() else output_candidate
        metadata = load_metadata_file(metadata_path)
        if metadata is not None:
            return metadata_path.resolve(), metadata

    return resolve_run_metadata(output_path, run_selector)


def resolve_status_log(output_path, metadata):
    snakemake_log_path = None
    fallback_logs = discover_snakemake_fallback_logs(output_path)
    if metadata:
        configured_log = metadata.get("snakemake_log")
        if configured_log:
            snakemake_log_path = Path(configured_log)
        elif metadata.get("run_id"):
            candidate = build_snakemake_log_path(output_path, metadata["run_id"])
            if candidate.exists():
                snakemake_log_path = candidate
    if snakemake_log_path is None and fallback_logs:
        snakemake_log_path = fallback_logs[0]
    return snakemake_log_path


def parse_wildcards(text):
    wildcards = {}
    for item in str(text or "").split(","):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        wildcards[key.strip()] = value.strip()
    return wildcards


def parse_rule_order(modules):
    workflow_dir = PACKAGE_DIR / "workflow"
    order = []
    seen = set()

    for module in modules or []:
        for rel_path in WORKFLOW_RULE_FILES.get(module, []):
            path = workflow_dir / rel_path
            try:
                text = path.read_text(encoding="utf-8")
            except OSError:
                continue
            for match in re.finditer(r"^\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", text, re.MULTILINE):
                name = match.group(1)
                if name not in seen:
                    order.append(name)
                    seen.add(name)

    return {name: index for index, name in enumerate(order)}


def _new_job_block(segment, rule_type, rule):
    return {
        "segment": segment,
        "rule_type": rule_type,
        "rule": rule,
        "jobid": None,
        "wildcards_raw": "",
        "wildcards": {},
        "input": "",
        "output": "",
        "started_order": None,
    }


def parse_snakemake_status(log_path, metadata=None):
    summary = {
        "rules": [],
        "samples": [],
        "overall": {
            "completed": 0,
            "total": 0,
            "started": 0,
            "failed": 0,
            "pending": 0,
            "percent": None,
        },
        "has_log": False,
    }
    if log_path is None or not Path(log_path).exists():
        return summary

    segment = 0
    segment_has_activity = False
    segment_planned = defaultdict(Counter)
    job_by_key = {}
    block_by_jobid = {}
    completed_keys = set()
    completed_without_block = Counter()
    failed_rules = Counter()
    started_order = 0
    current_block = None
    in_job_stats = False
    seen_job_stat_rows = False

    with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                if in_job_stats and seen_job_stat_rows:
                    in_job_stats = False
                    seen_job_stat_rows = False
                continue

            if stripped.startswith("Building DAG of jobs"):
                if segment_has_activity:
                    segment += 1
                    block_by_jobid = {}
                    current_block = None
                    segment_has_activity = False
                continue

            if stripped == "Job stats:":
                in_job_stats = True
                seen_job_stat_rows = False
                segment_has_activity = True
                continue

            if in_job_stats:
                if stripped.lower().startswith("job ") or set(stripped) <= {"-", " "}:
                    continue
                match = re.match(r"^([A-Za-z_][A-Za-z0-9_.-]*)\s+(\d+)\s*$", stripped)
                if match:
                    rule_name = match.group(1)
                    count = int(match.group(2))
                    if rule_name != "total":
                        segment_planned[segment][rule_name] = max(
                            segment_planned[segment][rule_name],
                            count,
                        )
                    seen_job_stat_rows = True
                    continue
                if seen_job_stat_rows:
                    in_job_stats = False
                    seen_job_stat_rows = False

            block_match = re.match(r"^(localrule|localcheckpoint|rule|checkpoint)\s+(.+?):\s*$", stripped)
            if block_match:
                segment_has_activity = True
                current_block = _new_job_block(segment, block_match.group(1), block_match.group(2))
                started_order += 1
                current_block["started_order"] = started_order
                continue

            if current_block is not None:
                jobid_match = re.match(r"^jobid:\s*(\d+)\s*$", stripped)
                if jobid_match:
                    current_block["jobid"] = jobid_match.group(1)
                    key = (segment, current_block["jobid"])
                    job_by_key[key] = current_block
                    block_by_jobid[current_block["jobid"]] = current_block
                    continue
                for field in ("wildcards", "input", "output"):
                    if stripped.startswith(f"{field}:"):
                        value = stripped.split(":", 1)[1].strip()
                        if field == "wildcards":
                            current_block["wildcards_raw"] = value
                            current_block["wildcards"] = parse_wildcards(value)
                        else:
                            current_block[field] = value
                        break

            finish_match = re.search(r"Finished jobid:\s*(\d+)(?:\s+\(Rule:\s*([^)]+)\))?", stripped)
            if finish_match:
                jobid = finish_match.group(1)
                rule_name = finish_match.group(2)
                key = (segment, jobid)
                completed_keys.add(key)
                if key not in job_by_key and rule_name:
                    completed_without_block[rule_name] += 1
                continue

            failed_match = re.match(r"(?:Error|RuleException|MissingOutputException|MissingInputException).*?\brule\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", stripped)
            if failed_match:
                failed_rules[failed_match.group(1)] += 1

    planned_counts = Counter()
    for counts in segment_planned.values():
        planned_counts.update(counts)

    started_counts = Counter()
    completed_counts = Counter(completed_without_block)
    samples_by_rule = defaultdict(set)
    assemblies_by_rule = defaultdict(set)
    jobs_by_rule = defaultdict(list)

    for key, job in job_by_key.items():
        rule = job["rule"]
        started_counts[rule] += 1
        jobs_by_rule[rule].append((key, job))
        wildcards = job.get("wildcards") or {}
        if wildcards.get("sample"):
            samples_by_rule[rule].add(wildcards["sample"])
        if wildcards.get("assembly"):
            assemblies_by_rule[rule].add(wildcards["assembly"])
        if key in completed_keys:
            completed_counts[rule] += 1

    all_rule_names = set(planned_counts) | set(started_counts) | set(completed_counts) | set(failed_rules)
    metadata_status = str((metadata or {}).get("status", "")).strip().lower()
    rules = []
    for rule in all_rule_names:
        total = max(planned_counts[rule], started_counts[rule], completed_counts[rule])
        completed = min(completed_counts[rule], total) if total else completed_counts[rule]
        if metadata_status == "success" and total:
            completed = total
        started = max(started_counts[rule], completed)
        failed = failed_rules[rule]
        if failed:
            status = "failed"
        elif total and completed >= total:
            status = "done"
        elif started > completed:
            status = "running"
        elif completed > 0:
            status = "partial"
        elif total:
            status = "pending"
        else:
            status = "observed"
        rules.append(
            {
                "rule": rule,
                "status": status,
                "completed": completed,
                "total": total,
                "started": started,
                "failed": failed,
                "samples": sorted(samples_by_rule[rule]),
                "assemblies": sorted(assemblies_by_rule[rule]),
            }
        )

    real_rules = [rule for rule in rules if rule["rule"] != "all"]
    overall_total = sum(rule["total"] for rule in real_rules)
    overall_completed = sum(min(rule["completed"], rule["total"]) for rule in real_rules)
    overall_started = sum(rule["started"] for rule in real_rules)
    overall_failed = sum(rule["failed"] for rule in real_rules)
    percent = int(round((overall_completed / overall_total) * 100)) if overall_total else None
    summary["overall"] = {
        "completed": overall_completed,
        "total": overall_total,
        "started": overall_started,
        "failed": overall_failed,
        "pending": max(0, overall_total - overall_completed),
        "percent": percent,
    }
    summary["rules"] = rules
    summary["jobs_by_rule"] = jobs_by_rule
    summary["completed_keys"] = completed_keys
    summary["has_log"] = True
    return summary


def read_json(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle) or {}
    except (OSError, json.JSONDecodeError):
        return {}


def load_sample_context(output_path):
    data_dir = Path(output_path) / "data"
    sample_names = set()
    for name in (
        "sample_to_reads1.json",
        "preprocessed_to_reads1.json",
        "transcriptome_to_reads1.json",
    ):
        sample_names.update(read_json(data_dir / name).keys())

    assembly_to_samples = read_json(data_dir / "assembly_to_samples.json")
    clean_assembly_map = {}
    for assembly, samples in assembly_to_samples.items():
        if isinstance(samples, list):
            clean_assembly_map[str(assembly)] = [str(sample) for sample in samples]
            sample_names.update(str(sample) for sample in samples)

    return sorted(sample_names), clean_assembly_map


def job_samples(job, sample_names, assembly_to_samples):
    wildcards = job.get("wildcards") or {}
    sample = wildcards.get("sample")
    if sample:
        return [sample]

    assembly = wildcards.get("assembly")
    if assembly:
        if assembly in assembly_to_samples:
            return assembly_to_samples[assembly]
        if assembly in sample_names:
            return [assembly]

    return []


def build_sample_rows(status_summary, sample_names, assembly_to_samples, rule_order, show_complete=False):
    jobs_by_rule = status_summary.get("jobs_by_rule") or {}
    completed_keys = status_summary.get("completed_keys") or set()
    sample_rule_stats = defaultdict(lambda: defaultdict(lambda: {"total": 0, "completed": 0, "failed": 0}))
    sample_rules_seen = defaultdict(set)

    failed_by_rule = {rule["rule"]: rule["failed"] for rule in status_summary.get("rules", [])}
    display_rules = {
        rule["rule"]
        for rule in status_summary.get("rules", [])
        if show_complete or not is_helper_rule(rule["rule"])
    }

    for rule, jobs in jobs_by_rule.items():
        if rule not in display_rules:
            continue
        for key, job in jobs:
            for sample in job_samples(job, sample_names, assembly_to_samples):
                stats = sample_rule_stats[sample][rule]
                stats["total"] += 1
                if key in completed_keys:
                    stats["completed"] += 1
                if failed_by_rule.get(rule):
                    stats["failed"] += 1
                sample_rules_seen[rule].add(sample)

    for rule, seen_samples in sample_rules_seen.items():
        if len(seen_samples) >= len(sample_names):
            continue
        planned = next((item["total"] for item in status_summary.get("rules", []) if item["rule"] == rule), 0)
        if planned >= len(sample_names):
            for sample in sample_names:
                sample_rule_stats[sample][rule]

    rows = []
    all_samples = sorted(set(sample_names) | set(sample_rule_stats))
    for sample in all_samples:
        rule_stats = sample_rule_stats.get(sample, {})
        ordered_rules = sorted(
            rule_stats,
            key=lambda name: (rule_order.get(name, 10_000), name),
        )
        total = len(ordered_rules)
        completed = 0
        failed = False
        current = None
        for rule in ordered_rules:
            stats = rule_stats[rule]
            rule_done = stats["total"] > 0 and stats["completed"] >= stats["total"]
            if rule_done:
                completed += 1
                continue
            if stats["failed"]:
                failed = True
            if current is None:
                current = rule
        if current is None and total:
            current = "complete"
        rows.append(
            {
                "sample": sample,
                "status": "failed" if failed else ("done" if total and completed >= total else "running" if completed else "pending"),
                "completed": completed,
                "total": total,
                "current": current or "not observed",
            }
        )
    return rows


def is_helper_rule(rule_name):
    return rule_name in HELPER_RULES


def progress_bar(completed, total, width=28):
    if not total:
        return "[" + "?" * width + "]"
    filled = int(round((completed / total) * width))
    filled = max(0, min(width, filled))
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def format_rule_detail(rule):
    details = []
    if rule.get("samples"):
        samples = rule["samples"]
        sample_text = ",".join(samples[:3])
        if len(samples) > 3:
            sample_text += f",+{len(samples) - 3}"
        details.append(f"samples={sample_text}")
    if rule.get("assemblies"):
        assemblies = rule["assemblies"]
        assembly_text = ",".join(assemblies[:3])
        if len(assemblies) > 3:
            assembly_text += f",+{len(assemblies) - 3}"
        details.append(f"assemblies={assembly_text}")
    if rule.get("failed"):
        details.append(f"failures={rule['failed']}")
    return "; ".join(details)


def print_run_header(output_path, metadata_path, metadata, log_path):
    section("DRAKKAR STATUS")
    print(f"Output directory: {output_path}")
    print(f"Run ID: {metadata.get('run_id', metadata_path.stem.removeprefix('drakkar_'))}")
    print(f"Command: {metadata.get('command', 'unknown')}")
    modules = metadata.get("modules") or []
    print(f"Modules: {', '.join(modules) if modules else 'unknown'}")
    print(f"Status: {metadata.get('status', 'unknown')}")
    if metadata.get("current_workflow"):
        print(f"Current workflow: {metadata['current_workflow']}")
    print(f"Started: {metadata.get('started_at', metadata.get('timestamp', 'unknown'))}")
    if metadata.get("finished_at"):
        print(f"Finished: {metadata['finished_at']}")
    print(f"Metadata file: {metadata_path}")
    print(f"Snakemake log: {log_path if log_path else 'not found'}")


def print_overall(status_summary):
    overall = status_summary["overall"]
    completed = overall["completed"]
    total = overall["total"]
    percent = overall["percent"]
    percent_text = f"{percent}%" if percent is not None else "unknown"
    section("OVERALL PROGRESS")
    print(f"Rule jobs: {progress_bar(completed, total)} {completed}/{total or '?'} ({percent_text})")
    print(
        "Started: "
        f"{overall['started']} | Pending: {overall['pending']} | Failed rule events: {overall['failed']}"
    )


def print_rule_rows(rules, rule_order, show_complete=False):
    display_rules = [
        rule for rule in rules
        if rule["rule"] != "all" and (show_complete or not is_helper_rule(rule["rule"]))
    ]
    display_rules.sort(key=lambda item: (rule_order.get(item["rule"], 10_000), item["rule"]))
    section("RULE STATUS")
    if not display_rules:
        print("No rule progress found in the Snakemake log.")
        return
    print(f"{'Rule':<34} {'Status':<9} {'Done':>9}  Details")
    for rule in display_rules:
        done = f"{rule['completed']}/{rule['total'] or '?'}"
        print(f"{rule['rule']:<34} {rule['status']:<9} {done:>9}  {format_rule_detail(rule)}")


def print_sample_rows(samples):
    section("SAMPLE STATUS")
    if not samples:
        print("No sample-specific rule progress found in the Snakemake log.")
        return
    print(f"{'Sample':<28} {'Status':<9} {'Stage':>9}  Current")
    for sample in samples:
        done = f"{sample['completed']}/{sample['total'] or '?'}"
        print(
            f"{sample['sample']:<28} {sample['status']:<9} "
            f"{done:>9}  {sample['current']}"
        )


def run_status(
    target=None,
    output_dir=None,
    run_id=None,
    show_complete=False,
    view="both",
):
    output_path, run_selector = resolve_status_target(target, output_dir, run_id)
    if not output_path.exists():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path is not a directory: {output_path}")
        return 1

    metadata_path, metadata = resolve_status_metadata(output_path, run_selector)
    if metadata is None:
        selector_text = f" matching {run_selector}" if run_selector else ""
        print(f"{ERROR}ERROR:{RESET} No Drakkar workflow run metadata{selector_text} found in {output_path}.")
        return 1

    log_path = resolve_status_log(output_path, metadata)
    modules = metadata.get("modules") or ([metadata.get("command")] if metadata.get("command") else [])
    rule_order = parse_rule_order(modules)
    status_summary = parse_snakemake_status(log_path, metadata=metadata)
    sample_names, assembly_to_samples = load_sample_context(output_path)
    sample_rows = build_sample_rows(
        status_summary,
        sample_names,
        assembly_to_samples,
        rule_order,
        show_complete=show_complete,
    )

    print_run_header(output_path, metadata_path, metadata, log_path)
    if not status_summary["has_log"]:
        print(f"{INFO}INFO:{RESET} No Snakemake log file found for this run.")
        return 0

    print_overall(status_summary)
    if view in {"both", "rules"}:
        print_rule_rows(status_summary["rules"], rule_order, show_complete=show_complete)
    if view in {"both", "samples"}:
        print_sample_rows(sample_rows)
    return 0
