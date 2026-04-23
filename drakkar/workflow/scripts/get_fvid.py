#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def parse_header(header):
    entry, _, description = header.partition(" ")
    entry = entry.strip()
    description = description.strip()

    vfc_match = re.search(r"\b(VFC\d+)\b", description)
    vfc = vfc_match.group(1) if vfc_match else ""

    vf_type = ""
    type_match = re.search(r"\[([^\[\]]+)\]\s*$", description)
    if type_match and type_match.group(1) != vfc:
        vf_type = type_match.group(1).strip()

    vf = description
    if not vf:
        vf = entry

    return {
        "entry": entry,
        "vf": vf,
        "vfc": vfc,
        "vf_type": vf_type,
    }


def iter_fasta_headers(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def main():
    parser = argparse.ArgumentParser(description="Create a VFDB mapping table from FASTA headers.")
    parser.add_argument("fasta", help="Path to VFDB FASTA file")
    parser.add_argument("output", help="Path to output TSV file")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("entry\tvf\tvfc\tvf_type\n")
        for header in iter_fasta_headers(args.fasta):
            record = parse_header(header)
            if record["entry"] in seen:
                continue
            seen.add(record["entry"])
            handle.write(
                "{entry}\t{vf}\t{vfc}\t{vf_type}\n".format(
                    entry=record["entry"],
                    vf=record["vf"],
                    vfc=record["vfc"],
                    vf_type=record["vf_type"],
                )
            )


if __name__ == "__main__":
    main()
