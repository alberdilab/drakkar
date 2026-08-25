"""Write a FASTA holding only the genes that no sequence-homology tool annotated.

Used to gate the Foldseek/ProstT5 structural-annotation step on the
homology-orphan genes: structure prediction is expensive, so we only fold the
proteins that KEGG/Pfam/CAZy (whichever are enabled) failed to hit. Prodigal
protein headers, HMMER tblout query ids, and dbCAN target names use the same
``{contig}_{orf}`` string, so genes can be matched directly by id.
"""

import argparse
import csv
from pathlib import Path

from Bio import SearchIO, SeqIO


DEFAULT_EVALUE_THRESHOLD = 1e-10
DBCAN_REQUIRED_COLUMNS = {"HMM Name", "Target Name", "i-Evalue", "Coverage"}


def has_content(path):
    if not path:
        return False
    path_obj = Path(path)
    return path_obj.is_file() and path_obj.stat().st_size > 0


def dbcan_gene_ids(path):
    """Return accepted target ids when path is a dbCAN HMM result table."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not DBCAN_REQUIRED_COLUMNS.issubset(reader.fieldnames):
            return None
        return {row["Target Name"] for row in reader if row.get("Target Name")}


def annotated_gene_ids(annotation_files, evalue_threshold):
    """Collect genes accepted by HMMER or dbCAN's native HMM filters."""
    annotated = set()
    for path in annotation_files:
        if not has_content(path):
            continue

        accepted_by_dbcan = dbcan_gene_ids(path)
        if accepted_by_dbcan is not None:
            annotated.update(accepted_by_dbcan)
            continue

        with open(path) as handle:
            for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
                for hit in queryresult.hits:
                    evalue = getattr(hit, "evalue", None)
                    if evalue is None or evalue <= evalue_threshold:
                        annotated.add(queryresult.id)
                        break
    return annotated


def main():
    parser = argparse.ArgumentParser(description="Extract genes without a sequence-homology annotation.")
    parser.add_argument("-faa", required=True, help="Prodigal protein FASTA for the MAG")
    parser.add_argument("-o", required=True, help="Output FASTA with the unannotated genes only")
    parser.add_argument(
        "-annotations", "-hmmer", dest="annotation_files", nargs="*", default=[],
        help="HMMER tblout or dbCAN HMM result files used to define 'annotated'",
    )
    parser.add_argument(
        "-evalue", "--evalue", type=float, default=DEFAULT_EVALUE_THRESHOLD,
        help="Maximum e-value for a hit to count as annotated. Default: 1e-10.",
    )
    args = parser.parse_args()

    annotated = annotated_gene_ids(args.annotation_files, args.evalue)

    output = Path(args.o)
    output.parent.mkdir(parents=True, exist_ok=True)

    if not has_content(args.faa):
        output.write_text("")
        return

    with open(output, "w") as out_handle:
        kept = (record for record in SeqIO.parse(args.faa, "fasta") if record.id not in annotated)
        SeqIO.write(kept, out_handle, "fasta")


if __name__ == "__main__":
    main()
