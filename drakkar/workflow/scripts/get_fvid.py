#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


VFC_PATTERN = re.compile(r"\b(VFC\d+)\b")
BRACKET_PATTERN = re.compile(r"\[([^\[\]]+)\]")
CLASSIFICATION_SEPARATOR_PATTERN = re.compile(r"\s+[-\u2013\u2014]\s+")
VFC_CATEGORY_PATTERN = re.compile(r"(?P<vf_type>.+?)\s*\((?P<vfc>VFC\d+)\)\s*$")
VFDB_MAPPING_SCHEMA = "drakkar-vfdb-v2"


def parse_vf_classification(description):
    """Return the VFDB category and VFC id embedded in a header description.

    Current VFDB Set B headers end with an organism block. The virulence-factor
    classification is the bracketed block containing the VFC identifier, for
    example ``[Phospholipase C (VF0470) - Exotoxin (VFC0235)]``.
    """
    for section in BRACKET_PATTERN.findall(description):
        if not VFC_PATTERN.search(section):
            continue

        # The category is the final component of the VFC-bearing section. This
        # also supports headers whose section is simply ``[Adherence (VFC0001)]``.
        category = CLASSIFICATION_SEPARATOR_PATTERN.split(section)[-1].strip()
        category_match = VFC_CATEGORY_PATTERN.fullmatch(category)
        if category_match:
            return category_match.group("vfc"), category_match.group("vf_type").strip()

    # Keep the VFC id when a non-standard header cannot safely yield a category.
    vfc_match = VFC_PATTERN.search(description)
    return (vfc_match.group(1) if vfc_match else ""), ""


def parse_header(header):
    entry, _, description = header.partition(" ")
    entry = entry.strip()
    description = description.strip()

    vfc, vf_type = parse_vf_classification(description)

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


def write_mapping(fasta_path, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("entry\tvf\tvfc\tvf_type\tmapping_schema\n")
        for header in iter_fasta_headers(fasta_path):
            record = parse_header(header)
            if record["entry"] in seen:
                continue
            seen.add(record["entry"])
            handle.write(
                "{entry}\t{vf}\t{vfc}\t{vf_type}\t{mapping_schema}\n".format(
                    entry=record["entry"],
                    vf=record["vf"],
                    vfc=record["vfc"],
                    vf_type=record["vf_type"],
                    mapping_schema=VFDB_MAPPING_SCHEMA,
                )
            )


def main():
    parser = argparse.ArgumentParser(description="Create a VFDB mapping table from FASTA headers.")
    parser.add_argument("fasta", help="Path to VFDB FASTA file")
    parser.add_argument("output", help="Path to output TSV file")
    args = parser.parse_args()
    write_mapping(args.fasta, args.output)


if __name__ == "__main__":
    main()
