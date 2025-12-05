#!/usr/bin/env python3
import csv
import argparse
from pathlib import Path
from collections import defaultdict, Counter


def parse_int(val):
    try:
        return int(str(val).replace(",", ""))
    except (TypeError, ValueError):
        return None


def load_cgc(file_path: Path):
    with file_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}

def extract_function_categories(label: str):
    if not label:
        return []
    categories = []
    for part in str(label).split("+"):
        head = part.split("|", 1)[0].strip()
        if head:
            categories.append(head)
    # Deduplicate preserving order
    return list(dict.fromkeys(categories))


def summarize_cgcs(rows):
    grouped = defaultdict(list)
    for row in rows:
        cgc_id = row.get("CGC#") or row.get("CGC") or row.get("cgc") or ""
        if not cgc_id:
            continue
        grouped[cgc_id].append(row)

    summaries = []
    for cgc_id, genes in grouped.items():
        contigs = [g.get("Contig ID") or g.get("contig") or "" for g in genes if (g.get("Contig ID") or g.get("contig"))]
        contig = contigs[0] if contigs else ""

        starts = [parse_int(g.get("Gene Start")) for g in genes if parse_int(g.get("Gene Start")) is not None]
        ends = [parse_int(g.get("Gene Stop")) for g in genes if parse_int(g.get("Gene Stop")) is not None]
        start = min(starts) if starts else ""
        end = max(ends) if ends else ""

        gene_types = [g.get("Gene Type") or "" for g in genes if g.get("Gene Type")]
        type_label = Counter(gene_types).most_common(1)[0][0] if gene_types else cgc_id

        annotations = []
        for g in genes:
            ann = g.get("Gene Annotation") or g.get("Gene annotation") or ""
            if not ann:
                continue
            cats = extract_function_categories(ann) or [ann]
            annotations.extend(cats)

        freq = Counter(annotations)
        annotated = [f"{ann} [{n}]" for ann, n in freq.items()]

        summaries.append({
            "contig": contig,
            "type": type_label,
            "start": start,
            "end": end,
            "gene_count": len(genes),
            "gene_functions": "; ".join(annotated),
        })

    return summaries


def write_summary(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["contig", "type", "start", "end", "gene_count", "gene_functions"],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize dbCAN cgc_standard_out.tsv into a cluster table.")
    parser.add_argument("-i", "--input", required=True, help="Path to cgc_standard_out.tsv")
    parser.add_argument("-o", "--output", required=True, help="Path to write summary TSV")
    return parser.parse_args()


def main():
    args = parse_args()
    rows = list(load_cgc(Path(args.input)))
    summary = summarize_cgcs(rows)
    write_summary(summary, Path(args.output))
    print(f"Wrote {len(summary)} clusters to {args.output}")


if __name__ == "__main__":
    main()
