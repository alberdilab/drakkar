#!/usr/bin/env python3
import csv
import argparse
from pathlib import Path


SUMMARY_COLUMNS = [
    "contig",
    "start",
    "end",
    "type",
    "gene_count",
    "substrate",
    "gene_functions",
    "pul_id",
    "source",
]


def read_table(path: Path, source: str, has_pul: bool = False):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({
                "contig": row.get("contig", ""),
                "start": row.get("start", ""),
                "end": row.get("end", ""),
                "type": row.get("type", ""),
                "gene_count": row.get("gene_count", ""),
                "substrate": row.get("substrate", ""),
                "gene_functions": row.get("gene_functions", ""),
                "pul_id": row.get("pul_id", "") if has_pul else "",
                "source": source,
            })
    return rows


def write_table(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Merge dbCAN, geNomad, and antiSMASH cluster summaries.")
    parser.add_argument("-dbcan", required=True, help="Path to dbCAN summary TSV")
    parser.add_argument("-genomad", required=True, help="Path to geNomad summary TSV")
    parser.add_argument("-antismash", required=True, help="Path to antiSMASH summary TSV")
    parser.add_argument("-o", "--output", required=True, help="Path to write merged TSV")
    return parser.parse_args()


def main():
    args = parse_args()
    merged = []
    merged.extend(read_table(Path(args.dbcan), "dbcan", has_pul=True))
    merged.extend(read_table(Path(args.genomad), "genomad"))
    merged.extend(read_table(Path(args.antismash), "antismash"))
    write_table(merged, Path(args.output))
    print(f"Wrote {len(merged)} merged clusters to {args.output}")


if __name__ == "__main__":
    main()
