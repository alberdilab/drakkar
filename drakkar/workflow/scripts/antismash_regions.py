#!/usr/bin/env python3
import re
import json
import sys
import csv
import argparse
from pathlib import Path
from html import unescape
from collections import Counter

def load_regions_json(js_path: Path):
    text = js_path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"=\s*(\[\s*[\s\S]*?\])\s*;?\s*$", text, re.MULTILINE)
    if not m:
        m = re.search(r"recordData\s*=\s*(\[\s*[\s\S]*?\])\s*;?\s*$", text, re.MULTILINE)
    if not m:
        raise ValueError("Could not locate JSON array in regions.js")
    json_text = m.group(1)

    try:
        return json.loads(json_text)
    except json.JSONDecodeError:
        relaxed = re.sub(r"(?<!\\)'", '"', json_text)
        relaxed = re.sub(r",(\s*[}\]])", r"\1", relaxed)
        return json.loads(relaxed)

def infer_region_type(region: dict) -> str | None:
    for k in ("type", "product", "products", "region_type", "most_significant"):
        val = region.get(k)
        if isinstance(val, list) and val:
            return str(val[0])
        if isinstance(val, str) and val.strip():
            return val.strip()

    orfs = region.get("orfs") or []
    hits = []
    for orf in orfs:
        desc = unescape(orf.get("description", "") or "")
        m = re.search(r"biosynthetic\s*\(rule-based-clusters\)\s*([A-Za-z0-9_-]+)", desc)
        if m:
            hits.append(m.group(1))
    if hits:
        return Counter(hits).most_common(1)[0][0]
    return None

def extract_orf_function(orf: dict) -> str | None:
    prod = (orf.get("product") or "").strip()
    if prod:
        return prod

    desc = unescape(orf.get("description", "") or "")
    if not desc:
        return None

    m = re.search(r"smcogs\)\s*SMCOG\d+:\s*([^<\n]+)", desc, flags=re.IGNORECASE)
    if m:
        return m.group(1).strip()

    m = re.search(r"\(rule-based-clusters\)\s*([A-Za-z0-9_.-]+)", desc, flags=re.IGNORECASE)
    if m:
        return m.group(1).strip()

    tail = re.findall(r">([^><\n]{2,})<\s*/div>", desc)
    if tail:
        guess = tail[-1].strip()
        guess = re.sub(r".*?:\s*", "", guess)
        if 2 <= len(guess) <= 60:
            return guess

    return "hypothetical protein"

def extract_function_categories(label: str):
    if not label:
        return []
    categories = []
    for part in str(label).split("+"):
        head = part.split("|", 1)[0].strip()
        if head:
            categories.append(head)
    return list(dict.fromkeys(categories))

def extract_gene_type(orf: dict) -> str:
    gene_functions = orf.get("gene_functions") or []
    labels = []
    for gf in gene_functions:
        if isinstance(gf, dict):
            label = gf.get("gene_function") or gf.get("function") or gf.get("category")
            if label:
                labels.append(str(label))
        elif isinstance(gf, str):
            labels.append(gf)
    if labels:
        return ";".join(dict.fromkeys(labels))

    for key in ("gene_kind", "kind", "type", "role"):
        val = orf.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()

    return ""

def format_strand(strand) -> str:
    if isinstance(strand, str) and strand.strip() in ("+", "-"):
        return strand.strip()
    if isinstance(strand, (int, float)):
        return "+" if strand > 0 else "-"
    return ""

def resolve_bgc_code(region: dict, contig: str, idx: int) -> str:
    for key in ("region_number", "region", "region_id", "region_name"):
        val = region.get(key)
        if isinstance(val, (int, float)):
            return f"BGC{int(val)}"
        if isinstance(val, str) and val.strip():
            return val.strip()
    return f"{contig}_BGC{idx}"

def write_summary(rows, out_csv: Path):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
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
            ]
        )
        writer.writeheader()
        writer.writerows(rows)

def write_gene_table(rows, out_path: Path):
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

def main(in_path: str, summary_csv: str, gene_table: str):
    records = load_regions_json(Path(in_path))
    summary_rows = []
    gene_rows = []
    bgc_index = 1

    for rec in records:
        contig = rec.get("seq_id") or rec.get("contig") or rec.get("record") or ""
        for region in rec.get("regions", []):
            start = region.get("start")
            end = region.get("end")
            rtype = infer_region_type(region) or ""
            bgc_code = resolve_bgc_code(region, contig, bgc_index)
            bgc_index += 1

            orfs = region.get("orfs") or []
            gene_count = len(orfs)

            # Count gene function frequencies
            funcs = []
            for o in orfs:
                annot = extract_orf_function(o) or ""
                cats = extract_function_categories(annot) or [annot]
                funcs.extend(cats)
            freq = Counter(funcs)
            # Build annotated list: func[count]
            annotated = [f"{f} [{n}]" for f, n in freq.items()]

            summary_rows.append({
                "contig": contig,
                "start": start,
                "end": end,
                "type": "BGC",
                "gene_count": gene_count,
                "gene_functions": "; ".join(annotated),
                "substrate": rtype,
            })

            for orf in orfs:
                protein_id = (
                    orf.get("protein_id")
                    or orf.get("locus_tag")
                    or orf.get("id")
                    or orf.get("name")
                    or ""
                )
                gene_rows.append({
                    "bgc": bgc_code,
                    "gene_type": extract_gene_type(orf),
                    "contig": contig,
                    "protein_id": protein_id,
                    "start": orf.get("start"),
                    "end": orf.get("end"),
                    "strand": format_strand(orf.get("strand")),
                    "annotation": extract_orf_function(orf) or "",
                })

    write_summary(summary_rows, Path(summary_csv))
    write_gene_table(gene_rows, Path(gene_table))

    print(f"Wrote {len(summary_rows)} regions to {summary_csv}")
    print(f"Wrote {len(gene_rows)} genes to {gene_table}")

def parse_args():
    parser = argparse.ArgumentParser(description="Parse antiSMASH regions.js into summary and gene tables.")
    parser.add_argument("-i", "--input", dest="input", help="Path to regions.js")
    parser.add_argument("-s", "--summary", dest="summary", help="Output CSV summarizing clusters.")
    parser.add_argument("-g", "--genes", dest="genes", help="Output TSV with per-gene annotations.")
    parser.add_argument("positional", nargs="*", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.input and args.summary and args.genes:
        return args.input, args.summary, args.genes

    if len(args.positional) == 3 and not (args.input or args.summary or args.genes):
        return args.positional[0], args.positional[1], args.positional[2]

    if len(args.positional) == 2 and not (args.input or args.summary or args.genes):
        # Backwards compatibility: derive gene table name next to summary output
        in_path, summary_path = args.positional
        gene_path = str(Path(summary_path).with_suffix("")) + "_genes.tsv"
        return in_path, summary_path, gene_path

    parser.error("Provide input, summary output, and gene output paths.")

if __name__ == "__main__":
    in_path, summary, genes = parse_args()
    main(in_path, summary, genes)
