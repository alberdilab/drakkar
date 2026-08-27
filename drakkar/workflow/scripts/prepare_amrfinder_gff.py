#!/usr/bin/env python3
"""Add the protein ``Name`` attributes required by AMRFinderPlus.

Prodigal's GFF ``ID`` and protein FASTA identifier agree, but its native GFF
does not include ``Name``.  AMRFinderPlus combined mode uses ``Name`` to join
the GFF to the protein FASTA, so preserve the native record and add it.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def convert(input_path, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    converted = 0
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
            protein_id = values.get("Name") or values.get("ID")
            if not protein_id:
                raise ValueError(
                    f"Prodigal GFF CDS line {line_number} has neither ID nor Name"
                )
            if "Name" not in values:
                attributes.append(("Name", protein_id))
            fields[8] = ";".join(
                key if value is None else f"{key}={value}"
                for key, value in attributes
            )
            destination.write("\t".join(fields) + "\n")
            converted += 1
    if converted == 0:
        raise ValueError(f"Prodigal GFF contains no feature records: {input_path}")
    return converted


def main():
    parser = argparse.ArgumentParser(
        description="Make a Prodigal GFF joinable to its protein FASTA by AMRFinderPlus."
    )
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
