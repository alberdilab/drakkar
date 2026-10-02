#!/usr/bin/env python3
"""Project Drakkar's lossless gene table into gifter's community input schema."""

from __future__ import annotations

import argparse
import csv
import json
import lzma
import re
from pathlib import Path


OUTPUT_COLUMNS = ["genome_id", "gene_id", "namespace", "accession"]

# This is a schema mapping, not a list of markers used by gifter. Every native
# accession from these Drakkar sources is projected; gifter remains responsible
# for deciding which namespace/accession pairs have biological meaning.
SOURCE_NAMESPACES = {
    "kegg": "KO",
    "cazy": "CAZY",
    "pfam": "PFAM",
    "ncbifam": "NCBIFAM",
    "tigrfam": "TIGRFAM",
}


def _open_text(path, mode):
    """Open plain or xz-compressed text without relying on filename guessing elsewhere."""
    path = Path(path)
    if path.suffix == ".xz":
        return lzma.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def _split_ec_values(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        values = []
        for item in value:
            values.extend(_split_ec_values(item))
        return values
    text = str(value).strip()
    if not text:
        return []
    values = []
    for item in re.split(r"[\s,;]+", text):
        item = item.strip().removeprefix("EC:").strip("[]")
        if item:
            values.append(item)
    return values


def ec_accessions(details):
    """Return EC accessions explicitly retained in one source's structured details."""
    try:
        payload = json.loads(details or "{}")
    except (TypeError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid annotation details JSON: {details!r}") from error
    if not isinstance(payload, dict):
        raise ValueError("Annotation details must be a JSON object")

    values = []
    values.extend(_split_ec_values(payload.get("ec")))
    values.extend(_split_ec_values(payload.get("ec_numbers")))

    profile_metadata = payload.get("profile_metadata")
    if isinstance(profile_metadata, dict):
        values.extend(_split_ec_values(profile_metadata.get("ec_numbers")))

    associations = payload.get("ec_associations")
    if isinstance(associations, list):
        for association in associations:
            if isinstance(association, dict):
                values.extend(_split_ec_values(association.get("ec")))

    return values


def row_markers(row):
    source = str(row.get("source") or "").strip().lower()
    annotation_id = str(row.get("annotation_id") or "").strip()
    namespace = SOURCE_NAMESPACES.get(source)
    if namespace and annotation_id:
        yield namespace, annotation_id
    for accession in ec_accessions(row.get("details")):
        yield "EC", accession


def project_gifter_input(input_path, output_path):
    """Write one stable, unique marker row per genome/gene/namespace/accession."""
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with _open_text(input_path, "rt") as source, _open_text(output_path, "wt") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"mag", "gene", "source", "annotation_id", "details"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Gene annotation table is missing required column(s): "
                + ", ".join(sorted(missing))
            )

        writer = csv.DictWriter(destination, delimiter="\t", fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()

        current_gene = None
        markers = set()
        output_rows = 0

        def flush():
            nonlocal output_rows
            if current_gene is None:
                return
            genome_id, gene_id = current_gene
            for namespace, accession in sorted(markers):
                writer.writerow({
                    "genome_id": genome_id,
                    "gene_id": gene_id,
                    "namespace": namespace,
                    "accession": accession,
                })
                output_rows += 1

        for row in reader:
            key = (str(row["mag"]), str(row["gene"]))
            if current_gene is not None and key != current_gene:
                flush()
                markers.clear()
            current_gene = key
            markers.update(row_markers(row))
        flush()

    return output_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="Long-form gene annotation TSV(.xz)")
    parser.add_argument("-o", "--output", required=True, help="gifter input TSV(.xz)")
    args = parser.parse_args()
    project_gifter_input(args.input, args.output)


if __name__ == "__main__":
    main()
