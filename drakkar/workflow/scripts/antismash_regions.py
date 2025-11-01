#!/usr/bin/env python3
import re
import json
import sys
import csv
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

def main(in_path: str, out_csv: str):
    records = load_regions_json(Path(in_path))
    rows = []

    for rec in records:
        contig = rec.get("seq_id") or rec.get("contig") or rec.get("record") or ""
        for region in rec.get("regions", []):
            start = region.get("start")
            end = region.get("end")
            rtype = infer_region_type(region) or ""

            orfs = region.get("orfs") or []
            gene_count = len(orfs)

            # Count gene function frequencies
            funcs = [extract_orf_function(o) for o in orfs]
            freq = Counter(funcs)
            # Build annotated list: func[count]
            annotated = [f"{f} [{n}]" for f, n in freq.items()]

            rows.append({
                "contig": contig,
                "type": rtype,
                "start": start,
                "end": end,
                "gene_count": gene_count,
                "gene_functions": "; ".join(annotated)
            })

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["contig", "type", "start", "end", "gene_count", "gene_functions"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} regions to {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python antismash_regions.py /path/to/regions.js output.csv", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
