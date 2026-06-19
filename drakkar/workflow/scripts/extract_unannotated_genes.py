"""Write a FASTA holding only the genes that no sequence-homology tool annotated.

Used to gate the Foldseek/ProstT5 structural-annotation step on the
homology-orphan genes: structure prediction is expensive, so we only fold the
proteins that KEGG/Pfam/CAZy (whichever are enabled) failed to hit. Prodigal
protein headers are ``>{contig}_{orf}`` and HMMER tblout query ids use the same
string, so genes can be matched directly by id.
"""

import argparse
from pathlib import Path

from Bio import SearchIO, SeqIO


DEFAULT_EVALUE_THRESHOLD = 1e-10


def has_content(path):
    if not path:
        return False
    path_obj = Path(path)
    return path_obj.is_file() and path_obj.stat().st_size > 0


def annotated_gene_ids(hmmer_files, evalue_threshold):
    """Collect query ids with at least one hit at or below the e-value cutoff."""
    annotated = set()
    for path in hmmer_files:
        if not has_content(path):
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
    parser.add_argument("-hmmer", nargs="*", default=[], help="HMMER tblout files (KEGG/Pfam/CAZy) used to define 'annotated'")
    parser.add_argument(
        "-evalue", "--evalue", type=float, default=DEFAULT_EVALUE_THRESHOLD,
        help="Maximum e-value for a hit to count as annotated. Default: 1e-10.",
    )
    args = parser.parse_args()

    annotated = annotated_gene_ids(args.hmmer, args.evalue)

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
