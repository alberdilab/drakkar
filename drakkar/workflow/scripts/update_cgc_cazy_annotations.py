#!/usr/bin/env python
import argparse
import csv
import os
import re
import sys
from typing import Dict, List


def load_cazy_types(hmm_results_path: str) -> Dict[str, List[str]]:
    """Return mapping from target (gene) to a list of CAZyme type prefixes."""
    type_regex = re.compile(r"^([A-Za-z]+)")
    gene_to_types: Dict[str, List[str]] = {}

    with open(hmm_results_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row.get("Target Name")
            hmm_name = row.get("HMM Name", "")
            if not gene:
                continue
            match = type_regex.match(hmm_name)
            if not match:
                continue
            cazy_type = match.group(1)
            gene_to_types.setdefault(gene, [])
            if cazy_type not in gene_to_types[gene]:
                gene_to_types[gene].append(cazy_type)

    return gene_to_types


def parse_attributes(attr_field: str) -> Dict[str, str]:
    attrs = {}
    for part in attr_field.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
        else:
            attrs[part] = ""
    return attrs


def attrs_to_string(attrs: Dict[str, str]) -> str:
    return ";".join(f"{k}={v}" if v != "" else k for k, v in attrs.items())


def update_cgc_annotations(gff_path: str, gene_to_types: Dict[str, List[str]]) -> None:
    updated_lines = []
    for line in open(gff_path):
        if line.startswith("#") or not line.strip():
            updated_lines.append(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            updated_lines.append(line)
            continue
        attrs = parse_attributes(parts[8])
        protein_id = attrs.get("protein_id")
        cgc_annotation = attrs.get("CGC_annotation", "null")

        types = gene_to_types.get(protein_id)
        if not types:
            updated_lines.append(line)
            continue

        existing = [] if cgc_annotation in ("null", "None", "") else cgc_annotation.split("|")
        for t in types:
            if t not in existing:
                existing.append(t)
        attrs["CGC_annotation"] = "|".join(existing) if existing else "null"
        parts[8] = attrs_to_string(attrs)
        updated_lines.append("\t".join(parts) + "\n")

    with open(gff_path, "w") as handle:
        handle.writelines(updated_lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add CAZyme type information from dbCAN hmm results into cgc.gff CGC_annotation field."
    )
    parser.add_argument("--hmm_results", required=True, help="Path to dbCAN_hmm_results.tsv")
    parser.add_argument("--gff", required=True, help="Path to cgc.gff to update in place")
    args = parser.parse_args()

    if not os.path.exists(args.hmm_results):
        sys.exit(f"ERROR: HMM results file not found: {args.hmm_results}")
    if not os.path.exists(args.gff):
        sys.exit(f"ERROR: GFF file not found: {args.gff}")

    gene_to_types = load_cazy_types(args.hmm_results)
    update_cgc_annotations(args.gff, gene_to_types)


if __name__ == "__main__":
    main()
