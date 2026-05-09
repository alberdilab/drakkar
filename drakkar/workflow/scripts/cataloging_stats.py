import argparse
import json
import os
import re

import pandas as pd


COLUMNS = [
    "assembly",
    "samples",
    "coverage_samples",
    "assembly_contigs",
    "assembly_total_length",
    "assembly_largest_contig",
    "assembly_gc_percent",
    "assembly_N50",
    "assembly_N75",
    "assembly_L50",
    "assembly_L75",
    "mapped_reads",
    "total_reads",
    "mapping_rate_percent",
    "sample_mapping_rates",
    "metabat2_bins",
    "maxbin2_bins",
    "semibin2_bins",
    "comebin_bins",
    "final_bins",
    "high_quality_bins",
    "medium_quality_bins",
    "low_quality_bins",
    "bin_total_size",
    "bin_mean_size",
    "bin_mean_N50",
    "bin_total_contigs",
    "bin_mean_completeness",
    "bin_mean_contamination",
    "best_bin",
    "best_bin_score",
    "best_bin_completeness",
    "best_bin_contamination",
    "best_bin_size",
    "best_bin_N50",
]

BINNER_METHODS = [
    ("metabat2_bins", "metabat2"),
    ("maxbin2_bins", "maxbin2"),
    ("semibin2_bins", "semibin2"),
    ("comebin_bins", "comebin"),
]

BINNER_ALIASES = {
    "metabat2_bins": ("metabat2", "metabat"),
    "maxbin2_bins": ("maxbin2", "maxbin"),
    "semibin2_bins": ("semibin2", "semibin"),
    "comebin_bins": ("comebin",),
}


def sample_from_path(path, suffix):
    name = os.path.basename(path)
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def clean_metric_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def numeric_or_na(value):
    if value in (None, "", "NA"):
        return "NA"
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number):
        return "NA"
    if float(number).is_integer():
        return int(number)
    return float(number)


def assembly_from_quast_path(path):
    return os.path.basename(os.path.dirname(path))


def assembly_sample_from_flagstat_path(path):
    sample = sample_from_path(path, ".flagstat.txt")
    assembly = os.path.basename(os.path.dirname(path))
    return assembly, sample


def read_quast_report(path):
    assembly = assembly_from_quast_path(path)
    metrics = {}

    try:
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                metrics[clean_metric_name(fields[0])] = fields[1]
    except OSError as error:
        print(f"Error processing {path}: {error}")

    def get_metric(*names):
        for name in names:
            value = metrics.get(clean_metric_name(name))
            if value not in (None, ""):
                return numeric_or_na(value)
        return "NA"

    return assembly, {
        "assembly_contigs": get_metric("# contigs", "# contigs (>= 0 bp)"),
        "assembly_total_length": get_metric("Total length", "Total length (>= 0 bp)"),
        "assembly_largest_contig": get_metric("Largest contig"),
        "assembly_gc_percent": get_metric("GC (%)"),
        "assembly_N50": get_metric("N50"),
        "assembly_N75": get_metric("N75"),
        "assembly_L50": get_metric("L50"),
        "assembly_L75": get_metric("L75"),
    }


def parse_count_pair(line):
    match = re.match(r"\s*(\d+)\s+\+\s+(\d+)\s+", line)
    if not match:
        return None
    return int(match.group(1)) + int(match.group(2))


def read_flagstat(path):
    assembly, sample = assembly_sample_from_flagstat_path(path)
    total_reads = 0
    mapped_reads = 0

    try:
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                if " in total " in line:
                    total_reads = parse_count_pair(line) or 0
                elif " mapped (" in line and "mate mapped" not in line:
                    mapped_reads = parse_count_pair(line) or 0
                    break
    except OSError as error:
        print(f"Error processing {path}: {error}")

    mapping_rate = (mapped_reads / total_reads * 100) if total_reads else 0
    return assembly, sample, {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "mapping_rate": mapping_rate,
    }


def read_tsv_row_count(path):
    if not path:
        return 0
    try:
        with open(path, "r", encoding="utf-8") as handle:
            lines = [line for line in handle if line.strip()]
    except OSError as error:
        print(f"Error processing {path}: {error}")
        return 0

    if not lines:
        return 0
    first_fields = {field.strip().lower() for field in lines[0].rstrip("\n").split("\t")}
    has_header = bool(first_fields.intersection({"bin_id", "contig", "contig_name", "file_name"}))
    return max(0, len(lines) - 1) if has_header else len(lines)


def read_quality_report_count(path):
    try:
        table = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
        print(f"Error processing {path}: {error}")
        return 0, None

    method = infer_binner_method(path, table)
    return len(table), method


def infer_binner_method(path, table):
    basename = os.path.basename(path)
    searchable = [basename, os.path.splitext(basename)[0]]
    for column in ("origin", "original_name", "name", "bin_id"):
        if column in table.columns:
            searchable.extend(str(value) for value in table[column].dropna().head(20))

    text = clean_metric_name(" ".join(searchable))
    for method, aliases in BINNER_ALIASES.items():
        if any(clean_metric_name(alias) in text for alias in aliases):
            return method
    return None


def binette_report_paths(report_root, assembly):
    report_dir = os.path.join(report_root, assembly, "input_bins_quality_reports")
    if not os.path.isdir(report_dir):
        return []

    paths = []
    for root, _, files in os.walk(report_dir):
        for filename in files:
            if filename.startswith("."):
                continue
            if filename.lower().endswith((".tsv", ".txt", ".csv")):
                paths.append(os.path.join(root, filename))
    return sorted(paths)


def read_binette_input_bin_counts(report_root, assembly, method_inputs):
    counts = {method: 0 for method, _ in BINNER_METHODS}
    report_paths = binette_report_paths(report_root, assembly)
    if not report_paths:
        return counts

    unmatched_counts = []
    assigned_methods = set()
    for path in report_paths:
        count, method = read_quality_report_count(path)
        if method in counts:
            counts[method] = count
            assigned_methods.add(method)
        else:
            unmatched_counts.append(count)

    valid_methods = [
        method
        for method, _ in BINNER_METHODS
        if read_tsv_row_count(method_inputs.get(method, "")) > 0
    ]
    if not valid_methods:
        valid_methods = [method for method, _ in BINNER_METHODS]

    unassigned_methods = [method for method in valid_methods if method not in assigned_methods]
    for method, count in zip(unassigned_methods, unmatched_counts):
        counts[method] = count

    return counts


def assembly_from_final_tsv(path):
    return os.path.splitext(os.path.basename(path))[0]


def read_bin_summary(path):
    assembly = assembly_from_final_tsv(path)
    empty = {
        "final_bins": 0,
        "high_quality_bins": 0,
        "medium_quality_bins": 0,
        "low_quality_bins": 0,
        "bin_total_size": 0,
        "bin_mean_size": "NA",
        "bin_mean_N50": "NA",
        "bin_total_contigs": 0,
        "bin_mean_completeness": "NA",
        "bin_mean_contamination": "NA",
        "best_bin": "NA",
        "best_bin_score": "NA",
        "best_bin_completeness": "NA",
        "best_bin_contamination": "NA",
        "best_bin_size": "NA",
        "best_bin_N50": "NA",
    }

    try:
        table = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
        print(f"Error processing {path}: {error}")
        return assembly, empty

    if table.empty:
        return assembly, empty

    for column in ["completeness", "contamination", "score", "size", "N50", "contig_count"]:
        if column in table.columns:
            table[column] = pd.to_numeric(table[column], errors="coerce")

    completeness = table["completeness"] if "completeness" in table.columns else pd.Series(dtype=float)
    contamination = table["contamination"] if "contamination" in table.columns else pd.Series(dtype=float)

    high_quality = (completeness >= 90) & (contamination <= 5)
    medium_quality = (completeness >= 50) & (contamination <= 10) & ~high_quality

    if "score" in table.columns and table["score"].notna().any():
        best_idx = table["score"].idxmax()
    elif completeness.notna().any():
        best_idx = completeness.idxmax()
    else:
        best_idx = table.index[0]
    best = table.loc[best_idx]

    summary = empty.copy()
    summary.update(
        {
            "final_bins": len(table),
            "high_quality_bins": int(high_quality.sum()) if not completeness.empty else 0,
            "medium_quality_bins": int(medium_quality.sum()) if not completeness.empty else 0,
            "low_quality_bins": int((~high_quality & ~medium_quality).sum()) if not completeness.empty else len(table),
            "bin_total_size": int(table["size"].sum()) if "size" in table.columns else "NA",
            "bin_mean_size": numeric_or_na(table["size"].mean()) if "size" in table.columns else "NA",
            "bin_mean_N50": numeric_or_na(table["N50"].mean()) if "N50" in table.columns else "NA",
            "bin_total_contigs": int(table["contig_count"].sum()) if "contig_count" in table.columns else "NA",
            "bin_mean_completeness": numeric_or_na(completeness.mean()) if not completeness.empty else "NA",
            "bin_mean_contamination": numeric_or_na(contamination.mean()) if not contamination.empty else "NA",
            "best_bin": str(best["bin_id"]) if "bin_id" in table.columns else "NA",
            "best_bin_score": numeric_or_na(best.get("score", "NA")),
            "best_bin_completeness": numeric_or_na(best.get("completeness", "NA")),
            "best_bin_contamination": numeric_or_na(best.get("contamination", "NA")),
            "best_bin_size": numeric_or_na(best.get("size", "NA")),
            "best_bin_N50": numeric_or_na(best.get("N50", "NA")),
        }
    )
    return assembly, summary


def load_assembly_samples(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        print(f"Error processing {path}: {error}")
        return {}

    return {
        str(assembly): ",".join(map(str, samples if isinstance(samples, list) else [samples]))
        for assembly, samples in data.items()
    }


def build_summary(args):
    assembly_samples = load_assembly_samples(args.assembly_to_samples)
    assemblies = set(assembly_samples)

    quast_data = {}
    for path in args.quast:
        assembly, data = read_quast_report(path)
        quast_data[assembly] = data
        assemblies.add(assembly)

    flagstat_data = {}
    for path in args.flagstat:
        assembly, sample, data = read_flagstat(path)
        flagstat_data.setdefault(assembly, {})[sample] = data
        assemblies.add(assembly)

    method_counts = {assembly: {} for assembly in assemblies}
    method_inputs = {}
    for method, paths in [
        ("metabat2_bins", args.metabat2),
        ("maxbin2_bins", args.maxbin2),
        ("semibin2_bins", args.semibin2),
        ("comebin_bins", args.comebin),
    ]:
        for path in paths:
            assembly = assembly_from_final_tsv(path)
            method_counts.setdefault(assembly, {})[method] = read_tsv_row_count(path)
            method_inputs.setdefault(assembly, {})[method] = path
            assemblies.add(assembly)

    bin_data = {}
    for path in args.bins:
        assembly, data = read_bin_summary(path)
        bin_data[assembly] = data
        assemblies.add(assembly)

    if args.binette_report_root:
        for assembly in assemblies:
            method_counts[assembly] = read_binette_input_bin_counts(
                args.binette_report_root,
                assembly,
                method_inputs.get(assembly, {}),
            )

    rows = []
    for assembly in sorted(assemblies):
        row = {column: "NA" for column in COLUMNS}
        row["assembly"] = assembly
        row["samples"] = assembly_samples.get(assembly, "NA")
        row.update(quast_data.get(assembly, {}))
        row.update(method_counts.get(assembly, {}))
        row.update(bin_data.get(assembly, {}))

        sample_stats = flagstat_data.get(assembly, {})
        if sample_stats:
            row["coverage_samples"] = ",".join(sorted(sample_stats))
            row["total_reads"] = sum(item["total_reads"] for item in sample_stats.values())
            row["mapped_reads"] = sum(item["mapped_reads"] for item in sample_stats.values())
            row["mapping_rate_percent"] = (
                row["mapped_reads"] / row["total_reads"] * 100 if row["total_reads"] else 0
            )
            row["sample_mapping_rates"] = ";".join(
                f"{sample}:{sample_stats[sample]['mapping_rate']:.2f}"
                for sample in sorted(sample_stats)
            )
        else:
            row["coverage_samples"] = "NA"
            row["total_reads"] = 0
            row["mapped_reads"] = 0
            row["mapping_rate_percent"] = 0
            row["sample_mapping_rates"] = "NA"

        rows.append(row)

    return pd.DataFrame(rows, columns=COLUMNS)


def main():
    parser = argparse.ArgumentParser(
        description="Create a cataloging summary table from assembly, mapping, and binning outputs."
    )
    parser.add_argument("--assembly-to-samples", required=True, help="assembly_to_samples JSON file")
    parser.add_argument("--quast", nargs="*", default=[], help="QUAST report.tsv files")
    parser.add_argument("--flagstat", nargs="*", default=[], help="samtools flagstat files")
    parser.add_argument("--metabat2", nargs="*", default=[], help="MetaBAT2 contig-to-bin TSV files")
    parser.add_argument("--maxbin2", nargs="*", default=[], help="MaxBin2 contig-to-bin TSV files")
    parser.add_argument("--semibin2", nargs="*", default=[], help="SemiBin2 contig-to-bin TSV files")
    parser.add_argument("--comebin", nargs="*", default=[], help="COMEBIN contig-to-bin TSV files")
    parser.add_argument(
        "--binette-report-root",
        help="Root cataloging/binette directory containing per-assembly input_bins_quality_reports directories",
    )
    parser.add_argument("--bins", nargs="*", default=[], help="Final per-assembly bin metadata TSV files")
    parser.add_argument("-o", "--output", required=True, help="Output cataloging TSV file")

    args = parser.parse_args()
    summary = build_summary(args)
    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
