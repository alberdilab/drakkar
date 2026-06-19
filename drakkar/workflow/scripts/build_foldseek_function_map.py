"""Build the Foldseek function-mapping table (FOLDSEEK_MAP_DB).

Foldseek searches the AlphaFold/Swiss-Prot database, whose target ids are
UniProt accessions. To turn a structural hit into a functional annotation we
need a UniProt accession -> KO / EC / Pfam table. UniProt's Swiss-Prot flat
file (``uniprot_sprot.dat.gz``) already carries those cross-references, so this
script parses it once into the TSV that ``merge_gene_annotations.py`` consumes.

The output has the header ``accession<TAB>kegg<TAB>ec<TAB>pfam``. When an entry
lists several KO/EC/Pfam values, the first of each is kept so the column stays a
single value, matching the sequence-homology annotation columns.

Stage the input with, e.g.::

    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

Then::

    python build_foldseek_function_map.py -i uniprot_sprot.dat.gz -o swissprot_function.tsv
"""

import argparse
import gzip
import re
from contextlib import contextmanager
from pathlib import Path

from Bio import SwissProt


EC_PATTERN = re.compile(r"EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9n-]+)")


@contextmanager
def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt")
    else:
        handle = open(path, "r")
    try:
        yield handle
    finally:
        handle.close()


def first_cross_reference(record, database):
    for reference in record.cross_references:
        if reference and reference[0] == database and len(reference) > 1:
            # Strip trailing punctuation that survives when the DR line has no
            # trailing description field (e.g. "DR   KO; K00010;").
            value = reference[1].strip().rstrip(";.").strip()
            if database == "Pfam":
                value = value.split(".")[0]
            return value
    return ""


def first_ec_number(record):
    match = EC_PATTERN.search(record.description)
    return match.group(1) if match else ""


def build_map(input_path, output_path):
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with open_text(input_path) as handle, open(output, "w") as out_handle:
        out_handle.write("accession\tkegg\tec\tpfam\n")
        for record in SwissProt.parse(handle):
            if not record.accessions:
                continue
            accession = record.accessions[0]
            kegg = first_cross_reference(record, "KO")
            ec = first_ec_number(record)
            pfam = first_cross_reference(record, "Pfam")
            # Skip entries that carry none of the three functional annotations.
            if not (kegg or ec or pfam):
                continue
            out_handle.write(f"{accession}\t{kegg}\t{ec}\t{pfam}\n")
            written += 1
    return written


def main():
    parser = argparse.ArgumentParser(description="Build the UniProt accession -> KO/EC/Pfam map for Foldseek hits.")
    parser.add_argument("-i", "--input", required=True, help="Path to uniprot_sprot.dat or uniprot_sprot.dat.gz")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path (FOLDSEEK_MAP_DB)")
    args = parser.parse_args()

    written = build_map(args.input, args.output)
    print(f"INFO Wrote {written} accession mappings to {args.output}")


if __name__ == "__main__":
    main()
