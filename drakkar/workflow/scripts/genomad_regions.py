#!/usr/bin/env python3
import csv
import argparse
from pathlib import Path
from collections import defaultdict, Counter


def parse_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def parse_coords(coords: str):
    if not coords or coords == "NA":
        return "", ""
    if "-" in coords:
        parts = coords.split("-", 1)
        try:
            return int(parts[0]), int(parts[1])
        except ValueError:
            return "", ""
    return "", ""


def format_strand(val):
    if val in ("+", "-"):
        return val
    try:
        f = float(val)
        return "+" if f > 0 else "-"
    except (TypeError, ValueError):
        return ""


def extract_function_categories(label: str):
    if not label:
        return []
    categories = []
    for part in str(label).split("+"):
        head = part.split("|", 1)[0].strip()
        if head:
            categories.append(head)
    return list(dict.fromkeys(categories))


def load_summary(path: Path, min_marker_enrichment: float, min_virus_score: float):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            marker = parse_float(row.get("marker_enrichment"))
            score = parse_float(row.get("virus_score"))
            if marker is not None and marker < min_marker_enrichment:
                continue
            if score is not None and score < min_virus_score:
                continue
            rows.append(row)
    return rows


def load_genes(path: Path):
    genes = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genes.append(row)
    return genes


def base_seq_name(gene_id: str):
    # Remove trailing _<number> to get parent sequence
    if not gene_id:
        return ""
    parts = gene_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return gene_id


def summarize(summary_rows, gene_rows):
    genes_by_seq = defaultdict(list)
    for g in gene_rows:
        parent = base_seq_name(g.get("gene") or g.get("seq_name") or "")
        genes_by_seq[parent].append(g)

    summaries = []
    gene_table = []

    for row in summary_rows:
        seq = row.get("seq_name") or ""
        coords = row.get("coordinates")
        start, end = parse_coords(coords)

        topology = (row.get("topology") or "").strip()
        type_label = "prophae" if topology == "Provirus" else "Viral fragment"

        gcount = row.get("n_genes")
        try:
            gcount_int = int(gcount)
        except (TypeError, ValueError):
            gcount_int = len(genes_by_seq.get(seq, []))

        taxonomy = row.get("taxonomy") or ""
        substrate = taxonomy if taxonomy else "NA"

        funcs = []
        for g in genes_by_seq.get(seq, []):
            ann = g.get("annotation_description") or g.get("annotation_accessions") or g.get("taxname") or ""
            cats = extract_function_categories(ann) or [ann] if ann else []
            funcs.extend(cats)
        freq = Counter(funcs)
        annotated = [f"{f} [{n}]" for f, n in freq.items()]

        summaries.append({
            "contig": seq,
            "start": start,
            "end": end,
            "type": type_label,
            "gene_count": gcount_int,
            "substrate": substrate,
            "gene_functions": "; ".join(annotated),
        })

        for g in genes_by_seq.get(seq, []):
            gene_table.append({
                "bgc": seq,
                "gene_type": g.get("marker") or "",
                "contig": seq,
                "protein_id": g.get("gene") or "",
                "start": g.get("start"),
                "end": g.get("end"),
                "strand": format_strand(g.get("strand")),
                "annotation": g.get("annotation_description") or g.get("annotation_accessions") or g.get("taxname") or "",
            })

    return summaries, gene_table


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
                "substrate",
                "gene_functions",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_genes(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "bgc",
                "gene_type",
                "contig",
                "protein_id",
                "start",
                "end",
                "strand",
                "annotation",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize geNomad virus outputs into region and gene tables.")
    parser.add_argument("-s", "--summary", required=True, help="Path to *_virus_summary.tsv")
    parser.add_argument("-g", "--genes", required=True, help="Path to *_virus_genes.tsv")
    parser.add_argument("-o", "--output", required=True, help="Path to write summary TSV")
    parser.add_argument("-G", "--gene-output", required=True, help="Path to write gene TSV")
    parser.add_argument("--min-marker-enrichment", type=float, default=0.0, help="Minimum marker_enrichment to keep entry")
    parser.add_argument("--min-virus-score", type=float, default=0.0, help="Minimum virus_score to keep entry")
    return parser.parse_args()


def main():
    args = parse_args()
    summary_rows = load_summary(Path(args.summary), args.min_marker_enrichment, args.min_virus_score)
    gene_rows = load_genes(Path(args.genes))
    summaries, gene_table = summarize(summary_rows, gene_rows)
    write_summary(summaries, Path(args.output))
    write_genes(gene_table, Path(args.gene_output))
    print(f"Wrote {len(summaries)} regions to {args.output}")
    print(f"Wrote {len(gene_table)} genes to {args.gene_output}")


if __name__ == "__main__":
    main()
