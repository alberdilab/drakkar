import os
import json
import pandas as pd
import argparse
import glob
import gzip

def extract_fastp_data(json_files):
    """Extract total reads before and after filtering from Fastp JSON files."""
    data_dict = {}

    for json_file in json_files:
        try:
            with open(json_file, "r") as f:
                report = json.load(f)

            sample_name = os.path.basename(json_file).replace(".json", "")
            reads_raw = report["summary"]["before_filtering"]["total_reads"]
            reads_postfiltering = report["summary"]["after_filtering"]["total_reads"]
            bases_raw = report["summary"]["before_filtering"]["total_bases"]
            bases_postfiltering = report["summary"]["after_filtering"]["total_bases"]

            data_dict[sample_name] = {
                "reads_raw": reads_raw,
                "bases_raw": bases_raw,
                "reads_discarded": reads_raw - reads_postfiltering,
                "bases_discarded": bases_raw - bases_postfiltering
            }

        except (KeyError, json.JSONDecodeError) as e:
            print(f"Error processing {json_file}: {e}")

    return data_dict

def extract_fastq_data(fastq_files):
    """Extract the number of reads and bases from FASTQ files."""
    data_dict = {}

    for fastq_file in fastq_files:
        sample_name = os.path.basename(fastq_file).replace("_1.fq.gz", "")
        read_count = 0
        base_count = 0

        try:
            with gzip.open(fastq_file, "rt") as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:  # Sequence line
                        read_count += 1
                        base_count += len(line.strip())

            data_dict[sample_name] = {
                "reads_metagenomic": read_count,
                "bases_metagenomic": base_count * 2 # doubled for also considering the reverse reads
            }

        except Exception as e:
            print(f"Error processing {fastq_file}: {e}")

    return data_dict

def extract_text_data(text_files, suffix):
    """Extract numbers from sample-specific text files."""
    data_dict = {}

    for text_file in text_files:
        try:
            sample_name = os.path.basename(text_file).replace(suffix, "")
            with open(text_file, "r") as f:
                value = f.readline().strip()  # Read the first line (expected number)
                data_dict[sample_name] = int(value)

        except Exception as e:
            print(f"Error processing {text_file}: {e}")

    return data_dict

def main():
    parser = argparse.ArgumentParser(description="Extract sequencing statistics from Fastp JSON, FASTQ, and BAM files")
    parser.add_argument("-p", "--fastp", required=True, nargs='+', help="Space-separated fastp JSON files")
    parser.add_argument("-f", "--fastq", required=False, nargs='+', help="Space-separated fastq files")
    parser.add_argument("-m", "--metagenomic_bases", required=False, nargs='+', help="Space-separated metabases files")
    parser.add_argument("-M", "--metagenomic_reads", required=False, nargs='+', help="Space-separated metareads files")
    parser.add_argument("-g", "--genomic_bases", required=False, nargs='+', help="Space-separated hostbases files")
    parser.add_argument("-G", "--genomic_reads", required=False, nargs='+', help="Space-separated hostreads files")
    parser.add_argument("-o", "--output", required=True, help="Output filename (e.g., 'fastp_summary.tsv')")

    args = parser.parse_args()

    # Collect data
    fastp_data = extract_fastp_data(args.fastp)
    fastq_data = extract_fastq_data(args.fastq) if args.fastq else {}
    metagenomic_bases_data = extract_text_data(args.metagenomic_bases,".metabases") if args.metagenomic_bases else {}
    metagenomic_reads_data = extract_text_data(args.metagenomic_reads,".metareads") if args.metagenomic_reads else {}
    genomic_bases_data = extract_text_data(args.genomic_bases,".hostbases") if args.genomic_bases else {}
    genomic_reads_data = extract_text_data(args.genomic_reads,".hostreads") if args.genomic_reads else {}

    # Combine data
    all_samples = set(fastp_data.keys()) | set(fastq_data.keys()) | set(metagenomic_bases_data.keys()) | set(metagenomic_reads_data.keys()) | set(genomic_bases_data.keys()) | set(genomic_reads_data.keys())
    summary_list = []

    for sample in all_samples:
        summary_list.append({
            "sample": sample,
            "reads_raw": fastp_data.get(sample, {}).get("reads_raw", "NA"),
            "reads_discarded": fastp_data.get(sample, {}).get("reads_discarded", "NA"),
            "reads_host": genomic_reads_data.get(sample, "NA"),
            "reads_metagenomic": fastq_data.get(sample, {}).get("reads_metagenomic", "NA"),
            "reads_metagenomic": metagenomic_reads_data.get(sample, "NA"),
            "bases_raw": fastp_data.get(sample, {}).get("bases_raw", "NA"),
            "bases_discarded": fastp_data.get(sample, {}).get("bases_discarded", "NA"),
            "bases_host": genomic_bases_data.get(sample, "NA"),
            "bases_metagenomic": fastq_data.get(sample, {}).get("bases_metagenomic", "NA"),
            "bases_metagenomic": metagenomic_bases_data.get(sample, "NA")
        })

    df = pd.DataFrame(summary_list)
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
