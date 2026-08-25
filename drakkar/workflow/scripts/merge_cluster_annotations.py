#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path


CLUSTER_QC_SCHEMA = "drakkar-cluster-annotation-qc-v1"
SUMMARY_COLUMNS = [
    "mag",
    "cluster_id",
    "contig",
    "start",
    "end",
    "source",
    "method",
    "evidence",
    "type",
    "gene_count",
    "substrate",
    "gene_functions",
    "pul_id",
    "details",
]

SOURCE_METADATA = {
    "dbcan": ("run_dbcan_cgc", "carbohydrate_gene_cluster_prediction"),
    "genomad": ("genomad", "mobile_element_prediction"),
    "antismash": ("antismash", "biosynthetic_gene_cluster_prediction"),
    "defensefinder": ("defense_finder", "defense_system_prediction"),
}

SOURCE_ALIASES = {
    "mobile": "genomad",
    "defense": "defensefinder",
}

GENERIC_REQUIRED_COLUMNS = {
    "cluster_id",
    "contig",
    "start",
    "end",
    "type",
    "gene_count",
    "substrate",
    "gene_functions",
}

DEFENSE_REQUIRED_COLUMNS = {
    "sys_id",
    "sys_beg",
    "sys_end",
    "type",
    "genes_count",
}


def has_content(path):
    if not path:
        return False
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def compact_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def validate_columns(path, fieldnames, required):
    missing = required.difference(fieldnames or [])
    if missing:
        raise ValueError(
            f"Cluster annotation input {path} is missing required columns: "
            f"{', '.join(sorted(missing))}. Regenerate this source with Drakkar 2.0."
        )


def fallback_cluster_id(source, row):
    identity = compact_json(row)
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:12]
    contig = row.get("contig") or "unknown-contig"
    start = row.get("start") or "unknown-start"
    end = row.get("end") or "unknown-end"
    return f"{source}:{contig}:{start}-{end}:{digest}"


def generic_rows(path, mag, source):
    if not path:
        return []
    path = Path(path)
    if not has_content(path):
        return []
    method, evidence = SOURCE_METADATA[source]
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(path, reader.fieldnames, GENERIC_REQUIRED_COLUMNS)
        for native in reader:
            cluster_id = native.get("cluster_id") or fallback_cluster_id(source, native)
            rows.append({
                "mag": str(mag),
                "cluster_id": cluster_id,
                "contig": native.get("contig", ""),
                "start": native.get("start", ""),
                "end": native.get("end", ""),
                "source": source,
                "method": method,
                "evidence": evidence,
                "type": native.get("type", ""),
                "gene_count": native.get("gene_count", ""),
                "substrate": native.get("substrate", ""),
                "gene_functions": native.get("gene_functions", ""),
                "pul_id": native.get("pul_id", "") if source == "dbcan" else "",
                "details": compact_json(native),
            })
    return rows


def defense_rows(path, mag):
    if not path:
        return []
    path = Path(path)
    if not has_content(path):
        return []
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(path, reader.fieldnames, DEFENSE_REQUIRED_COLUMNS)
        for native in reader:
            cluster_id = native.get("sys_id") or fallback_cluster_id("defensefinder", native)
            contig = next(
                (
                    native.get(column, "")
                    for column in ("replicon", "contig", "sequence", "seq_id")
                    if native.get(column)
                ),
                "",
            )
            rows.append({
                "mag": str(mag),
                "cluster_id": cluster_id,
                "contig": contig,
                "start": native.get("sys_beg", ""),
                "end": native.get("sys_end", ""),
                "source": "defensefinder",
                "method": "defense_finder",
                "evidence": "defense_system_prediction",
                "type": native.get("type", ""),
                "gene_count": native.get("genes_count", ""),
                "substrate": native.get("subtype", ""),
                "gene_functions": native.get("name_of_profiles_in_sys", ""),
                "pul_id": "",
                "details": compact_json(native),
            })
    return rows


def normalize_sources(sources):
    if sources is None:
        return {"dbcan", "genomad", "antismash", "defensefinder"}
    if isinstance(sources, str):
        sources = sources.split(",")
    normalized = {
        SOURCE_ALIASES.get(str(source).strip().lower(), str(source).strip().lower())
        for source in sources
        if str(source).strip()
    }
    unknown = normalized.difference(SOURCE_METADATA)
    if unknown:
        raise ValueError(f"Unknown enabled cluster annotation sources: {', '.join(sorted(unknown))}")
    return normalized


def validate_cluster_ids(rows):
    seen = set()
    duplicates = set()
    for row in rows:
        key = (row["mag"], row["source"], row["cluster_id"])
        if key in seen:
            duplicates.add(":".join(key))
        seen.add(key)
    if duplicates:
        raise ValueError(
            "Cluster annotations contain duplicate (MAG, source, cluster_id) keys: "
            + ", ".join(sorted(duplicates)[:10])
        )


def write_table(rows, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def write_qc(rows_by_source, mag, out_path):
    records = []
    for source, rows in rows_by_source.items():
        records.append({
            "mag": str(mag),
            "level": "cluster",
            "source": source,
            "reported_records": len(rows),
            "retained_records": len(rows),
            "rejected_records": None,
            "unmapped_records": sum(not row.get("contig") for row in rows),
            "unique_entities": len({row["cluster_id"] for row in rows}),
            "filter_stage": "upstream_native",
        })
    payload = {
        "schema_version": CLUSTER_QC_SCHEMA,
        "mag": str(mag),
        "level": "cluster",
        "sources": records,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(compact_json(payload) + "\n", encoding="utf-8")


def merge_cluster_annotations(
    *,
    mag,
    output,
    sources=None,
    dbcan=None,
    genomad=None,
    antismash=None,
    defense=None,
    qc_output=None,
):
    if not mag:
        raise ValueError("MAG identity is required for the cluster annotation table")
    enabled = normalize_sources(sources)
    rows_by_source = {}
    if "dbcan" in enabled:
        rows_by_source["dbcan"] = generic_rows(dbcan, mag, "dbcan")
    if "genomad" in enabled:
        rows_by_source["genomad"] = generic_rows(genomad, mag, "genomad")
    if "antismash" in enabled:
        rows_by_source["antismash"] = generic_rows(antismash, mag, "antismash")
    if "defensefinder" in enabled:
        rows_by_source["defensefinder"] = defense_rows(defense, mag)

    merged = [row for rows in rows_by_source.values() for row in rows]
    validate_cluster_ids(merged)
    write_table(merged, output)
    if qc_output:
        write_qc(rows_by_source, mag, qc_output)
    return merged


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge validated dbCAN, geNomad, antiSMASH, and DefenseFinder clusters."
    )
    parser.add_argument("-dbcan", required=False, help="Path to dbCAN summary TSV")
    parser.add_argument("-genomad", required=False, help="Path to geNomad summary TSV")
    parser.add_argument("-antismash", required=False, help="Path to antiSMASH summary TSV")
    parser.add_argument("-defense", required=False, help="Path to DefenseFinder systems TSV")
    parser.add_argument("--mag", required=True, help="MAG identifier written to every row")
    parser.add_argument(
        "--sources",
        default=None,
        help="Comma-separated enabled sources; existing files from other sources are ignored.",
    )
    parser.add_argument("--qc-output", help="Path to write per-MAG annotation QC JSON")
    parser.add_argument("-o", "--output", required=True, help="Path to write merged TSV")
    return parser.parse_args()


def main():
    args = parse_args()
    rows = merge_cluster_annotations(
        mag=args.mag,
        output=args.output,
        sources=args.sources,
        dbcan=args.dbcan,
        genomad=args.genomad,
        antismash=args.antismash,
        defense=args.defense,
        qc_output=args.qc_output,
    )
    print(f"Wrote {len(rows)} merged clusters to {args.output}")


if __name__ == "__main__":
    main()
