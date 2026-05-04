#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_reads(path):
    line_count = 0
    total_bases = 0
    with open_maybe_gzip(path) as handle:
        for line_count, line in enumerate(handle, start=1):
            if line_count % 4 == 2:
                total_bases += len(line.strip())
    return line_count // 4, total_bases


def main():
    parser = argparse.ArgumentParser(description="Count total reads across one or more FASTQ files.")
    parser.add_argument("fastq", nargs="+", help="Input FASTQ or FASTQ.GZ files")
    parser.add_argument("-o", "--output", required=True, help="Output file containing the total read count")
    parser.add_argument(
        "-b",
        "--bases-output",
        required=False,
        help="Optional output file containing the total number of sequenced bases",
    )
    args = parser.parse_args()

    total_reads = 0
    total_bases = 0
    for path in args.fastq:
        reads, bases = count_reads(path)
        total_reads += reads
        total_bases += bases

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(f"{total_reads}\n", encoding="utf-8")

    if args.bases_output:
        bases_output_path = Path(args.bases_output)
        bases_output_path.parent.mkdir(parents=True, exist_ok=True)
        bases_output_path.write_text(f"{total_bases}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
