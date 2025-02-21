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
    parser.add_argument("-p", "--fastp", required=True, help="Glob pattern for Fastp JSON files (e.g., 'fastp_results/*.json')")
    parser.add_argument("-f", "--fastq", required=False, help="Glob pattern for FASTQ files (e.g., 'fastq_files/*_1.fq.gz')")
    parser.add_argument("-m", "--metagenomic_bases", required=False, help="Glob pattern for FASTQ files (e.g., 'fastq_files/*_1.fq.gz')")
    parser.add_argument("-M", "--metagenomic_reads", required=False, help="Glob pattern for FASTQ files (e.g., 'fastq_files/*_1.fq.gz')")
    parser.add_argument("-g", "--genomic_bases", required=False, help="Glob pattern for BAM files (e.g., 'bam_files/*.bam')")
    parser.add_argument("-G", "--genomic_reads", required=False, help="Glob pattern for BAM files (e.g., 'bam_files/*.bam')")
    parser.add_argument("-o", "--output", required=True, help="Output filename (e.g., 'fastp_summary.tsv')")

    args = parser.parse_args()

    # Collect data
    json_files = glob.glob(args.fastp)
    fastp_data = extract_fastp_data(json_files)
    if args.fastq:
        fastq_files = glob.glob(args.fastq)
        fastq_data = extract_fastq_data(fastq_files)
    else:
        fastq_data = {}
    if args.metagenomic_bases:
        metagenomic_bases_files = glob.glob(args.metagenomic_bases)
        metagenomic_bases_data = extract_text_data(metagenomic_bases_files)
    else:
        metagenomic_bases_data = {}
    if args.metagenomic_reads:
        metagenomic_reads_files = glob.glob(args.metagenomic_reads)
        metagenomic_reads_data = extract_text_data(metagenomic_reads_files)
    else:
        metagenomic_reads_data = {}
    if args.genomic_bases:
        genomic_bases_files = glob.glob(args.genomic_bases)
        genomic_bases_data = extract_text_data(genomic_bases_files)
    else:
        genomic_bases_data = {}
    if args.genomic_reads:
        genomic_reads_files = glob.glob(args.genomic_reads)
        genomic_reads_data = extract_text_data(genomic_reads_files)
    else:
        genomic_reads_data = {}

    # Combine data
    all_samples = set(fastp_data.keys()) | set(metareads_data.keys()) | set(metabases_data.keys()) | set(hostreads_data.keys()) | set(hostbases_data.keys())
    summary_list = []

    for sample in all_samples:
        summary_list.append({
            "sample": sample,
            "reads_raw": fastp_data.get(sample, {}).get("reads_raw", "NA"),
            "bases_raw": fastp_data.get(sample, {}).get("bases_raw", "NA"),
            "reads_discarded": fastp_data.get(sample, {}).get("reads_discarded", "NA"),
            "bases_discarded": fastp_data.get(sample, {}).get("bases_discarded", "NA"),
            "reads_metagenomic": fastq_data.get(sample, {}).get("reads_metagenomic", "NA"),
            "bases_metagenomic": fastq_data.get(sample, {}).get("bases_metagenomic", "NA"),
            "reads_metagenomic": metagenomic_reads_data.get(sample, "NA"),
            "bases_metagenomic": metagenomic_bases_data.get(sample, "NA"),
            "reads_host": genomic_reads_data.get(sample, "NA"),
            "bases_host": genomic_bases_data.get(sample, "NA")
        })

    df = pd.DataFrame(summary_list)
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
