import argparse
import csv
import os

OUTPUT_COLUMNS = ["contig", "bin", "assembly"]

FASTA_EXTENSIONS = (".fa", ".fna", ".fasta")


def bin_name(fasta_path):
    """Return the bin identifier of a FASTA file (its basename without extension)."""
    name = os.path.basename(fasta_path)
    for extension in FASTA_EXTENSIONS:
        if name.endswith(extension):
            return name[: -len(extension)]
    return name


def contig_names(fasta_path):
    """Yield the identifier of every sequence in a FASTA file."""
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                header = line[1:].strip()
                if header:
                    yield header.split()[0]


def build_table(fasta_files, assembly, output_file):
    """Write the contig-to-bin table of one assembly's final bins."""
    rows = 0
    with open(output_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for fasta_path in fasta_files:
            genome = bin_name(fasta_path)
            for contig in contig_names(fasta_path):
                writer.writerow([contig, genome, assembly])
                rows += 1
    print(f"✅ Contig-to-bin table with {rows} contigs saved to {output_file}")


def merge_tables(table_files, output_file):
    """Concatenate per-assembly contig-to-bin tables into a single CSV."""
    rows = 0
    with open(output_file, "w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for table_file in table_files:
            with open(table_file, newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if header is not None and header != OUTPUT_COLUMNS:
                    raise ValueError(
                        f"Unexpected columns in {table_file}: {header}"
                    )
                for row in reader:
                    writer.writerow(row)
                    rows += 1
    print(f"✅ Combined contig-to-bin table with {rows} contigs saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map the contigs of the final bins back to the assembly they were binned from."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="Bin FASTA files, or per-assembly contig-to-bin tables when --merge is used",
    )
    parser.add_argument("--assembly", help="Assembly the bins were generated from")
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Combine per-assembly tables instead of reading bin FASTA files",
    )
    parser.add_argument("-o", "--output", required=True, help="Output table path")

    args = parser.parse_args()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if args.merge:
        merge_tables(args.inputs, args.output)
    else:
        if not args.assembly:
            parser.error("--assembly is required unless --merge is used")
        build_table(args.inputs, args.assembly, args.output)
