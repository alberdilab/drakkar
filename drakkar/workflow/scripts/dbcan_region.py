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

def load_substrates(file_path: Path):
    """
    Returns mapping of CGC id (e.g., CGC4) to tuple (PULID, substrate).
    """
    mapping = {}
    with file_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            raw_cgc = row.get("#cgcid") or row.get("cgcid") or row.get("CGC") or ""
            if not raw_cgc:
                continue
            cgc_id = str(raw_cgc).split("|")[-1].strip()
            if not cgc_id:
                continue
            pul_id = (row.get("PULID") or "").strip()
            substrate = (row.get("dbCAN-PUL substrate") or row.get("substrate") or "").strip()
            mapping[cgc_id] = (pul_id, substrate)
    return mapping

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


def summarize_cgcs(rows, pul_map=None):
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
        type_label = "CGC"
        pul_id = ""
        substrate = ""
        if pul_map:
            # Match on CGC name (e.g., CGC4)
            hit = pul_map.get(cgc_id)
            if hit:
                type_label = "PUL"
                pul_id, substrate = hit

        if type_label != "PUL":
            type_label = "CGC"

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
            "start": start,
            "end": end,
            "type": type_label,
            "gene_count": len(genes),
            "gene_functions": "; ".join(annotated),
            "substrate": substrate,
            "pul_id": pul_id,
        })

    return summaries


def write_summary(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "contig",
                "start",
                "end",
                "type",
                "gene_count",
                "gene_functions",
                "substrate",
                "pul_id",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize dbCAN cgc_standard_out.tsv into a cluster table.")
    parser.add_argument("-i", "--input", required=True, help="Path to cgc_standard_out.tsv")
    parser.add_argument("-o", "--output", required=True, help="Path to write summary TSV")
    parser.add_argument("-p", "--puls", help="Path to substrate_prediction.tsv (optional)")
    return parser.parse_args()


def main():
    args = parse_args()
    rows = list(load_cgc(Path(args.input)))
    pul_map = load_substrates(Path(args.puls)) if args.puls else None
    summary = summarize_cgcs(rows, pul_map=pul_map)
    write_summary(summary, Path(args.output))
    print(f"Wrote {len(summary)} clusters to {args.output}")


if __name__ == "__main__":
    main()
