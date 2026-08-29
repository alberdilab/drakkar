#!/usr/bin/env python3
"""Add the protein identifiers required by AMRFinderPlus to a Prodigal GFF.

Prodigal names each protein in its FASTA output ``<seqid>_<gene ordinal>``,
but its GFF ``ID`` attribute is ``<sequence ordinal>_<gene ordinal>`` (e.g.
``1_1``) and it writes no ``Name``.  AMRFinderPlus combined mode joins the GFF
to the protein FASTA by that identifier, so rebuild it from the GFF sequence
column and write it to both ``ID`` and ``Name``.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def protein_ids(proteins_path):
    """Return the identifiers of a Prodigal protein FASTA, in file order."""
    identifiers = []
    with open(proteins_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split(None, 1)[0])
    return identifiers


def convert(input_path, output_path, proteins_path=None):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    converted = 0
    names = []
    with open(input_path, encoding="utf-8") as source, output_path.open(
        "w", encoding="utf-8"
    ) as destination:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Prodigal GFF line {line_number} has {len(fields)} fields instead of 9"
                )
            attributes = []
            values = {}
            for item in fields[8].split(";"):
                if not item:
                    continue
                key, separator, value = item.partition("=")
                attributes.append((key, value if separator else None))
                if separator:
                    values[key] = value
            gene_id = values.get("ID")
            if not gene_id:
                raise ValueError(
                    f"Prodigal GFF CDS line {line_number} has no ID attribute"
                )
            # Prodigal's ID is <sequence ordinal>_<gene ordinal>; the protein
            # FASTA uses the sequence name instead of its ordinal.
            protein_id = f"{fields[0]}_{gene_id.rsplit('_', 1)[-1]}"
            attributes = [
                (key, protein_id if key in ("ID", "Name") else value)
                for key, value in attributes
            ]
            if "Name" not in values:
                attributes.append(("Name", protein_id))
            fields[8] = ";".join(
                key if value is None else f"{key}={value}"
                for key, value in attributes
            )
            destination.write("\t".join(fields) + "\n")
            names.append(protein_id)
            converted += 1
    if converted == 0:
        raise ValueError(f"Prodigal GFF contains no feature records: {input_path}")
    if proteins_path is not None:
        missing = sorted(set(protein_ids(proteins_path)) - set(names))
        if missing:
            raise ValueError(
                f"{len(missing)} protein(s) of {proteins_path} are absent from "
                f"{input_path}, first: {missing[0]}"
            )
    return converted


def main():
    parser = argparse.ArgumentParser(
        description="Make a Prodigal GFF joinable to its protein FASTA by AMRFinderPlus."
    )
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument(
        "--proteins",
        help="Prodigal protein FASTA to check the rewritten identifiers against",
    )
    args = parser.parse_args()
    convert(args.input, args.output, args.proteins)


if __name__ == "__main__":
    main()
