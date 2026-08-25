#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import yaml


MANIFEST_SCHEMA = "drakkar-annotation-manifest-v1"
QC_COLUMNS = [
    "mag",
    "level",
    "source",
    "reported_records",
    "retained_records",
    "rejected_records",
    "unmapped_records",
    "unique_entities",
    "filter_stage",
]


def parse_assignments(values):
    parsed = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"Expected NAME=VALUE, received: {value}")
        key, item = value.split("=", 1)
        parsed[key] = item
    return parsed


def sha256sum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_version():
    package_init = Path(__file__).resolve().parents[2] / "__init__.py"
    match = re.search(
        r'^__version__\s*=\s*["\']([^"\']+)["\']',
        package_init.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    return match.group(1) if match else "unknown"


def locate_database_manifest(configured_path):
    path = Path(configured_path)
    candidates = [
        path / "database_versions.yaml",
        path.parent / "database_versions.yaml",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def describe_database(configured_path):
    path = Path(configured_path)
    record = {
        "configured_path": str(path),
        "exists": path.exists(),
    }
    version_manifest = locate_database_manifest(path)
    if version_manifest:
        installed = yaml.safe_load(version_manifest.read_text(encoding="utf-8")) or {}
        record.update({
            "installation_manifest": str(version_manifest),
            "requested_version": installed.get("requested_version"),
            "source_version": installed.get("source_version"),
            "sources": installed.get("sources", []),
            "files": installed.get("files", []),
        })
    else:
        record["installation_manifest"] = None
        record["release_hint"] = path.name
        record["files"] = []
    return record


def describe_environment(path):
    path = Path(path)
    record = {"path": str(path), "exists": path.is_file(), "dependencies": []}
    if path.is_file():
        environment = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        record["dependencies"] = environment.get("dependencies", [])
    return record


def describe_table(path):
    path = Path(path)
    return {
        "path": str(path),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256sum(path) if path.is_file() else None,
    }


def load_qc_records(paths):
    records = []
    for path in paths:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        records.extend(payload.get("sources", []))
    return sorted(records, key=lambda row: (row.get("mag", ""), row.get("level", ""), row.get("source", "")))


def write_qc_table(records, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=QC_COLUMNS)
        writer.writeheader()
        for record in records:
            writer.writerow({column: record.get(column, "") for column in QC_COLUMNS})


def write_annotation_report(
    *,
    qc_inputs,
    manifest_output,
    qc_output,
    enabled_sources,
    thresholds=None,
    databases=None,
    tools=None,
    environments=None,
    tables=None,
):
    qc_records = load_qc_records(qc_inputs)
    write_qc_table(qc_records, qc_output)

    manifest = {
        "schema_version": MANIFEST_SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "drakkar_version": package_version(),
        "enabled_sources": sorted(enabled_sources),
        "thresholds": thresholds or {},
        "configured_tools": tools or {},
        "environments": {
            name: describe_environment(path) for name, path in (environments or {}).items()
        },
        "databases": {
            name: describe_database(path) for name, path in (databases or {}).items()
        },
        "tables": {
            name: describe_table(path) for name, path in (tables or {}).items()
        },
        "qc_table": str(Path(qc_output)),
    }
    output = Path(manifest_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return manifest, qc_records


def parse_args():
    parser = argparse.ArgumentParser(description="Write annotation provenance and QC sidecars.")
    parser.add_argument("--qc-input", nargs="+", required=True)
    parser.add_argument("--manifest-output", required=True)
    parser.add_argument("--qc-output", required=True)
    parser.add_argument("--enabled-sources", required=True)
    parser.add_argument("--threshold", nargs="*", default=[])
    parser.add_argument("--database", nargs="*", default=[])
    parser.add_argument("--tool", nargs="*", default=[])
    parser.add_argument("--environment", nargs="*", default=[])
    parser.add_argument("--table", nargs="*", default=[])
    return parser.parse_args()


def main():
    args = parse_args()
    write_annotation_report(
        qc_inputs=args.qc_input,
        manifest_output=args.manifest_output,
        qc_output=args.qc_output,
        enabled_sources=[source for source in args.enabled_sources.split(",") if source],
        thresholds=parse_assignments(args.threshold),
        databases=parse_assignments(args.database),
        tools=parse_assignments(args.tool),
        environments=parse_assignments(args.environment),
        tables=parse_assignments(args.table),
    )


if __name__ == "__main__":
    main()
