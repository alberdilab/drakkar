#!/usr/bin/env python3
"""
profiling_genomes_stats.py

Script to aggregate per-sample mapped reads and mapped bases counts
from lists of files into a single TSV summary.
"""
import argparse
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate per-sample mapped reads and bases into one summary TSV"
    )
    parser.add_argument(
        '-r', '--mappedreads',
        nargs='+',
        required=True,
        help='List of .mappedreads files'
    )
    parser.add_argument(
        '-b', '--mappedbases',
        nargs='+',
        required=True,
        help='List of .mappedbases files'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output TSV file (columns: sample, reads_mapped, bases_mapped)'
    )
    return parser.parse_args()


def read_counts(files, suffix_name):
    counts = {}
    for filepath in files:
        path = Path(filepath)
        if not path.exists():
            print(f"Warning: file not found: {filepath}", file=sys.stderr)
            continue
        sample = path.stem
        try:
            val = path.read_text().strip()
        except Exception as e:
            print(f"Error reading {suffix_name} file '{filepath}': {e}", file=sys.stderr)
            continue
        counts[sample] = val
    return counts


def main():
    args = parse_args()

    # Read counts into dictionaries keyed by sample
    reads_dict = read_counts(args.mappedreads, 'reads')
    bases_dict = read_counts(args.mappedbases, 'bases')

    # Determine all sample names
    samples = sorted(set(reads_dict) | set(bases_dict))
    if not samples:
        sys.exit("No samples found in provided files.")

    # Write summary
    with open(args.output, 'w') as out_f:
        out_f.write("sample\treads_mapped\tbases_mapped\n")
        for sample in samples:
            reads = reads_dict.get(sample, 'NA')
            bases = bases_dict.get(sample, 'NA')
            if reads == 'NA':
                print(f"Warning: reads count missing for sample '{sample}'", file=sys.stderr)
            if bases == 'NA':
                print(f"Warning: bases count missing for sample '{sample}'", file=sys.stderr)
            out_f.write(f"{sample}\t{reads}\t{bases}\n")

if __name__ == '__main__':
    main()
