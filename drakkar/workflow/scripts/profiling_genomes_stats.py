#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate per-sample profiling statistics into one summary TSV."
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
        '-t', '--totalreads',
        nargs='+',
        required=True,
        help='List of .totalreads files'
    )
    parser.add_argument(
        '-B', '--totalbases',
        nargs='+',
        required=True,
        help='List of .totalbases files'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output TSV file'
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


def parse_int(value, label, sample):
    try:
        return int(str(value))
    except ValueError:
        print(f"Warning: invalid {label} value for sample '{sample}': {value}", file=sys.stderr)
        return None


def main():
    args = parse_args()

    # Read counts into dictionaries keyed by sample
    reads_dict = read_counts(args.mappedreads, 'reads')
    bases_dict = read_counts(args.mappedbases, 'bases')
    total_reads_dict = read_counts(args.totalreads, 'total reads')
    total_bases_dict = read_counts(args.totalbases, 'total bases')

    # Determine all sample names
    samples = sorted(set(reads_dict) | set(bases_dict) | set(total_reads_dict) | set(total_bases_dict))
    if not samples:
        sys.exit("No samples found in provided files.")

    # Write summary
    with open(args.output, 'w') as out_f:
        out_f.write("sample\tinput_reads\tinput_bases\treads_mapped\tbases_mapped\tmapping_percentage\n")
        for sample in samples:
            reads_raw = reads_dict.get(sample, 'NA')
            bases_raw = bases_dict.get(sample, 'NA')
            total_reads_raw = total_reads_dict.get(sample, 'NA')
            total_bases_raw = total_bases_dict.get(sample, 'NA')
            if reads_raw == 'NA':
                print(f"Warning: reads count missing for sample '{sample}'", file=sys.stderr)
            if bases_raw == 'NA':
                print(f"Warning: bases count missing for sample '{sample}'", file=sys.stderr)
            if total_reads_raw == 'NA':
                print(f"Warning: total read count missing for sample '{sample}'", file=sys.stderr)
            if total_bases_raw == 'NA':
                print(f"Warning: total base count missing for sample '{sample}'", file=sys.stderr)

            reads = parse_int(reads_raw, "mapped reads", sample) if reads_raw != 'NA' else None
            bases = parse_int(bases_raw, "mapped bases", sample) if bases_raw != 'NA' else None
            total_reads = parse_int(total_reads_raw, "total reads", sample) if total_reads_raw != 'NA' else None
            total_bases = parse_int(total_bases_raw, "total bases", sample) if total_bases_raw != 'NA' else None

            mapping_percentage = "NA"
            if reads is not None and total_reads not in (None, 0):
                mapping_percentage = f"{(reads / total_reads) * 100:.2f}"

            out_f.write(
                f"{sample}\t"
                f"{total_reads if total_reads is not None else 'NA'}\t"
                f"{total_bases if total_bases is not None else 'NA'}\t"
                f"{reads if reads is not None else 'NA'}\t"
                f"{bases if bases is not None else 'NA'}\t"
                f"{mapping_percentage}\n"
            )

if __name__ == '__main__':
    main()
