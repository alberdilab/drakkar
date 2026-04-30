import argparse
import gzip
import json
import os

import pandas as pd


COLUMNS = [
    "sample",
    "reads_pre_fastp",
    "bases_pre_fastp",
    "adapter_trimmed_reads",
    "adapter_trimmed_bases",
    "reads_post_fastp",
    "bases_post_fastp",
    "host_reads",
    "host_bases",
    "metagenomic_reads",
    "metagenomic_bases",
    "singlem_fraction",
    "nonpareil_C",
    "nonpareil_LR",
    "nonpareil_modelR",
    "nonpareil_LRstar",
    "nonpareil_diversity",
]


def sample_from_path(path, suffix):
    name = os.path.basename(path)
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def first_present(mapping, keys, default="NA"):
    for key in keys:
        if key in mapping:
            return mapping[key]
    return default


def extract_fastp_data(json_files):
    data = {}

    for json_file in json_files:
        try:
            with open(json_file, "r", encoding="utf-8") as handle:
                report = json.load(handle)

            sample = sample_from_path(json_file, ".json")
            before = report["summary"]["before_filtering"]
            after = report["summary"]["after_filtering"]
            adapters = report.get("adapter_cutting", {})

            data[sample] = {
                "reads_pre_fastp": before.get("total_reads", "NA"),
                "bases_pre_fastp": before.get("total_bases", "NA"),
                "adapter_trimmed_reads": adapters.get("adapter_trimmed_reads", "NA"),
                "adapter_trimmed_bases": adapters.get("adapter_trimmed_bases", "NA"),
                "reads_post_fastp": after.get("total_reads", "NA"),
                "bases_post_fastp": after.get("total_bases", "NA"),
            }
        except (KeyError, json.JSONDecodeError, OSError) as error:
            print(f"Error processing {json_file}: {error}")

    return data


def extract_fastq_data(fastq_files):
    data = {}

    for fastq_file in fastq_files:
        sample = sample_from_path(fastq_file, "_1.fq.gz")
        read_count = 0
        base_count = 0

        try:
            with gzip.open(fastq_file, "rt") as handle:
                for line_number, line in enumerate(handle):
                    if line_number % 4 == 1:
                        read_count += 1
                        base_count += len(line.strip())

            data[sample] = {
                "metagenomic_reads": read_count * 2,
                "metagenomic_bases": base_count * 2,
            }
        except OSError as error:
            print(f"Error processing {fastq_file}: {error}")

    return data


def extract_text_data(text_files, suffix):
    data = {}

    for text_file in text_files:
        sample = sample_from_path(text_file, suffix)
        try:
            with open(text_file, "r", encoding="utf-8") as handle:
                value = handle.readline().strip()
            data[sample] = int(value)
        except (OSError, ValueError) as error:
            print(f"Error processing {text_file}: {error}")

    return data


def extract_singlem_fraction(singlem_files):
    data = {}
    candidate_columns = [
        "read_fraction",
        "fraction",
        "microbial_fraction",
        "singlem_fraction",
        "estimated_fraction",
    ]

    for singlem_file in singlem_files:
        fallback_sample = sample_from_path(singlem_file, "_smf.tsv")
        try:
            table = pd.read_csv(singlem_file, sep="\t")
            if table.empty:
                data[fallback_sample] = "NA"
                continue

            sample = table["sample"].iloc[0] if "sample" in table.columns else fallback_sample
            fraction = "NA"
            for column in candidate_columns:
                if column in table.columns:
                    fraction = table[column].iloc[0]
                    break

            if fraction == "NA":
                numeric_columns = table.select_dtypes(include="number").columns
                if len(numeric_columns) > 0:
                    fraction = table[numeric_columns[0]].iloc[0]

            data[str(sample)] = fraction
        except (OSError, pd.errors.ParserError, ValueError) as error:
            print(f"Error processing {singlem_file}: {error}")

    return data


def extract_nonpareil_data(nonpareil_files):
    data = {}

    for nonpareil_file in nonpareil_files:
        fallback_sample = sample_from_path(nonpareil_file, "_np.tsv")
        try:
            table = pd.read_csv(nonpareil_file, sep="\t")
            if table.empty:
                data[fallback_sample] = {}
                continue

            row = table.iloc[0]
            sample = row["sample"] if "sample" in table.columns else fallback_sample
            data[str(sample)] = {
                "nonpareil_C": first_present(row, ["C"]),
                "nonpareil_LR": first_present(row, ["LR"]),
                "nonpareil_modelR": first_present(row, ["modelR", "modelRt"]),
                "nonpareil_LRstar": first_present(row, ["LRstar"]),
                "nonpareil_diversity": first_present(row, ["diversity"]),
            }
        except (OSError, pd.errors.ParserError, ValueError) as error:
            print(f"Error processing {nonpareil_file}: {error}")

    return data


def build_summary(args):
    fastp_data = extract_fastp_data(args.fastp)
    fastq_data = extract_fastq_data(args.fastq) if args.fastq else {}
    metagenomic_bases_data = (
        extract_text_data(args.metagenomic_bases, ".metabases")
        if args.metagenomic_bases
        else {}
    )
    metagenomic_reads_data = (
        extract_text_data(args.metagenomic_reads, ".metareads")
        if args.metagenomic_reads
        else {}
    )
    host_bases_data = (
        extract_text_data(args.host_bases, ".hostbases") if args.host_bases else {}
    )
    host_reads_data = (
        extract_text_data(args.host_reads, ".hostreads") if args.host_reads else {}
    )
    singlem_data = (
        extract_singlem_fraction(args.singlem_fraction) if args.singlem_fraction else {}
    )
    nonpareil_data = extract_nonpareil_data(args.nonpareil) if args.nonpareil else {}

    all_samples = (
        set(fastp_data)
        | set(fastq_data)
        | set(metagenomic_bases_data)
        | set(metagenomic_reads_data)
        | set(host_bases_data)
        | set(host_reads_data)
        | set(singlem_data)
        | set(nonpareil_data)
    )

    rows = []
    for sample in sorted(all_samples):
        row = {column: "NA" for column in COLUMNS}
        row["sample"] = sample
        row.update(fastp_data.get(sample, {}))
        row.update(fastq_data.get(sample, {}))

        if sample in metagenomic_reads_data:
            row["metagenomic_reads"] = metagenomic_reads_data[sample]
        if sample in metagenomic_bases_data:
            row["metagenomic_bases"] = metagenomic_bases_data[sample]
        if sample in host_reads_data:
            row["host_reads"] = host_reads_data[sample]
        if sample in host_bases_data:
            row["host_bases"] = host_bases_data[sample]
        if sample in singlem_data:
            row["singlem_fraction"] = singlem_data[sample]
        row.update(nonpareil_data.get(sample, {}))

        rows.append(row)

    return pd.DataFrame(rows, columns=COLUMNS)


def main():
    parser = argparse.ArgumentParser(
        description="Extract preprocessing statistics from fastp, FASTQ, SingleM, and Nonpareil outputs"
    )
    parser.add_argument("-p", "--fastp", required=True, nargs="+", help="fastp JSON files")
    parser.add_argument("-f", "--fastq", required=False, nargs="+", help="R1 FASTQ files")
    parser.add_argument("-m", "--metagenomic_bases", required=False, nargs="+", help="metagenomic base count files")
    parser.add_argument("-M", "--metagenomic_reads", required=False, nargs="+", help="metagenomic read count files")
    parser.add_argument("-g", "--host_bases", "--genomic_bases", dest="host_bases", required=False, nargs="+", help="host base count files")
    parser.add_argument("-G", "--host_reads", "--genomic_reads", dest="host_reads", required=False, nargs="+", help="host read count files")
    parser.add_argument("-s", "--singlem_fraction", required=False, nargs="+", help="SingleM microbial fraction TSV files")
    parser.add_argument("-n", "--nonpareil", required=False, nargs="+", help="Nonpareil summary TSV files")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    args = parser.parse_args()
    summary = build_summary(args)
    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
