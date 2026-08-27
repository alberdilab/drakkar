#!/usr/bin/env python3
"""Aggregate per-assembly AMR digests and record reproducible provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
from datetime import datetime, timezone
from pathlib import Path

import yaml

from amr_digest import (
    DIGEST_COLUMNS,
    DRUG_COLUMNS,
    HIT_COLUMNS,
    LOCUS_COLUMNS,
    MOBILITY_COLUMNS,
    REGION_COLUMNS,
)


TABLES = {
    "amr_hits.tsv.xz": ("amr_hits.tsv", HIT_COLUMNS),
    "amr_loci.tsv.xz": ("amr_loci.tsv", LOCUS_COLUMNS),
    "amr_drug_classes.tsv.xz": ("amr_drug_classes.tsv", DRUG_COLUMNS),
    "mobility_regions.tsv.xz": ("mobility_regions.tsv", REGION_COLUMNS),
    "amr_mobility.tsv.xz": ("amr_mobility.tsv", MOBILITY_COLUMNS),
}

SUMMARY_COLUMNS = [
    "assembly_id", "assembly_type", "organism", "input_path", "input_sha256",
    "input_size_bytes", "contig_count", "total_length", "amrfinder_hits",
    "rgi_hits", "mobility_regions", "amr_loci", "multi_tool_loci",
    "mobile_loci",
]

QC_COLUMNS = [
    "assembly_id", "amrfinder_hits", "amrfinder_hits_without_coordinates",
    "rgi_hits", "rgi_hits_without_coordinates", "mobility_regions", "amr_loci",
    "multi_tool_loci", "mobility_links", "mobile_loci",
]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def database_record(configured):
    path = Path(configured)
    return {"configured": str(path), "release": path.name}


def read_rows(path, expected_columns):
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_columns:
            raise ValueError(
                f"Unexpected columns in {path}: {reader.fieldnames}; expected {expected_columns}"
            )
        return list(reader)


def write_rows(path, columns, rows, compressed=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = lzma.open if compressed else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=columns, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def concatenate_tables(output_path, input_paths, columns):
    """Stream per-assembly tables into one XZ table with bounded memory."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with lzma.open(output_path, "wt", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(
            output, delimiter="\t", fieldnames=columns, extrasaction="ignore"
        )
        writer.writeheader()
        for input_path in input_paths:
            with Path(input_path).open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != columns:
                    raise ValueError(
                        f"Unexpected columns in {input_path}: {reader.fieldnames}; "
                        f"expected {columns}"
                    )
                for row in reader:
                    writer.writerow(row)
                    count += 1
    return count


def source_count(qc, source, key):
    for record in qc.get("sources", []):
        if record.get("source") == source:
            return int(record.get(key, 0))
    return 0


def aggregate(args):
    manifest_path = Path(args.input_manifest)
    assemblies = json.loads(manifest_path.read_text(encoding="utf-8"))
    per_assembly = Path(args.per_assembly_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    qc_rows = []
    for assembly_id in sorted(assemblies):
        assembly_dir = per_assembly / assembly_id
        qc = json.loads((assembly_dir / "qc.json").read_text(encoding="utf-8"))
        amrfinder_hits = source_count(qc, "amrfinderplus", "retained_hits")
        amrfinder_no_coords = source_count(
            qc, "amrfinderplus", "hits_without_coordinates"
        )
        rgi_hits = source_count(qc, "rgi", "retained_hits")
        rgi_no_coords = source_count(qc, "rgi", "hits_without_coordinates")
        regions = source_count(qc, "genomad", "retained_regions")
        metadata = assemblies[assembly_id]
        qc_row = {
            "assembly_id": assembly_id,
            "amrfinder_hits": amrfinder_hits,
            "amrfinder_hits_without_coordinates": amrfinder_no_coords,
            "rgi_hits": rgi_hits,
            "rgi_hits_without_coordinates": rgi_no_coords,
            "mobility_regions": regions,
            "amr_loci": int(qc.get("loci", 0)),
            "multi_tool_loci": int(qc.get("multi_tool_loci", 0)),
            "mobility_links": int(qc.get("mobility_links", 0)),
            "mobile_loci": int(qc.get("mobile_loci", 0)),
        }
        qc_rows.append(qc_row)
        summaries.append({
            "assembly_id": assembly_id,
            "assembly_type": metadata.get("assembly_type", ""),
            "organism": metadata.get("organism", ""),
            "input_path": metadata.get("path", ""),
            "input_sha256": metadata.get("sha256", ""),
            "input_size_bytes": metadata.get("size_bytes", ""),
            "contig_count": metadata.get("contig_count", ""),
            "total_length": metadata.get("total_length", ""),
            "amrfinder_hits": amrfinder_hits,
            "rgi_hits": rgi_hits,
            "mobility_regions": regions,
            "amr_loci": qc_row["amr_loci"],
            "multi_tool_loci": qc_row["multi_tool_loci"],
            "mobile_loci": qc_row["mobile_loci"],
        })

    table_counts = {}
    for output_name, (input_name, columns) in TABLES.items():
        table_counts[output_name] = concatenate_tables(
            output_dir / output_name,
            [per_assembly / assembly_id / input_name for assembly_id in sorted(assemblies)],
            columns,
        )
    write_rows(output_dir / "assembly_summary.tsv", SUMMARY_COLUMNS, summaries)
    write_rows(output_dir / "amr_qc.tsv", QC_COLUMNS, qc_rows)

    described_outputs = {}
    for name in [*TABLES, "assembly_summary.tsv", "amr_qc.tsv"]:
        path = output_dir / name
        row_count = (
            table_counts[name] if name in table_counts
            else len(summaries) if name == "assembly_summary.tsv"
            else len(qc_rows)
        )
        described_outputs[name] = {
            "path": str(path.resolve()),
            "rows": row_count,
            "sha256": sha256(path),
        }

    provenance = {
        "schema_version": "drakkar-amr-manifest-v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "drakkar_version": args.drakkar_version,
        "assemblies": assemblies,
        "parameters": {
            "minimum_locus_overlap": args.minimum_overlap,
            "rgi_alignment_tool": args.rgi_alignment_tool,
            "rgi_include_loose": args.rgi_include_loose,
            "rgi_include_nudge": args.rgi_include_nudge,
            "genomad_preset": args.genomad_preset,
            "genomad_splits": args.genomad_splits,
        },
        "software": {
            "amrfinderplus": args.amrfinder_version,
            "prodigal": args.prodigal_version,
            "rgi": args.rgi_version,
            "genomad": args.genomad_version,
        },
        "databases": {
            "AMRFINDER_DB": database_record(args.amrfinder_db),
            "CARD_DB": database_record(args.card_db),
            "GENOMAD_DB": database_record(args.genomad_db),
        },
        "outputs": described_outputs,
    }
    (output_dir / "manifest.yaml").write_text(
        yaml.safe_dump(provenance, sort_keys=False), encoding="utf-8"
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--per-assembly-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--minimum-overlap", type=float, required=True)
    parser.add_argument("--rgi-alignment-tool", required=True)
    parser.add_argument("--rgi-include-loose", action="store_true")
    parser.add_argument("--rgi-include-nudge", action="store_true")
    parser.add_argument("--genomad-preset", required=True)
    parser.add_argument("--genomad-splits", type=int, required=True)
    parser.add_argument("--amrfinder-db", required=True)
    parser.add_argument("--card-db", required=True)
    parser.add_argument("--genomad-db", required=True)
    parser.add_argument("--drakkar-version", default="unknown")
    parser.add_argument("--amrfinder-version", default="4.2.7")
    parser.add_argument("--prodigal-version", default="2.6.3")
    parser.add_argument("--rgi-version", default="6.0.8")
    parser.add_argument("--genomad-version", default="1.11.0")
    return parser.parse_args()


if __name__ == "__main__":
    aggregate(parse_args())
