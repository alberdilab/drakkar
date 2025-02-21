import os
import json
import pandas as pd
import argparse
import glob
import gzip
import pysam

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
                "bases_metagenomic": base_count
            }

        except Exception as e:
            print(f"Error processing {fastq_file}: {e}")

    return data_dict

def extract_bam_data(bam_files):
    """Extract the number of reads and bases from BAM files."""
    data_dict = {}

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).replace(".bam", "")
        read_count = 0
        base_count = 0

        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam:
                    read_count += 1
                    base_count += read.query_length if read.query_length else 0

            data_dict[sample_name] = {
                "reads_genomic": read_count,
                "bases_genomic": base_count
            }

        except Exception as e:
            print(f"Error processing {bam_file}: {e}")

    return data_dict

def main():
    parser = argparse.ArgumentParser(description="Extract sequencing statistics from Fastp JSON, FASTQ, and BAM files")
    parser.add_argument("-f", "--fastp", required=True, help="Glob pattern for Fastp JSON files (e.g., 'fastp_results/*.json')")
    parser.add_argument("-m", "--metagenomic", required=True, help="Glob pattern for FASTQ files (e.g., 'fastq_files/*_1.fq.gz')")
    parser.add_argument("-b", "--bam", required=True, help="Glob pattern for BAM files (e.g., 'bam_files/*.bam')")
    parser.add_argument("-o", "--output", required=True, help="Output filename (e.g., 'fastp_summary.tsv')")

    args = parser.parse_args()

    # Collect data
    json_files = glob.glob(args.fastp)
    fastq_files = glob.glob(args.metagenomic)
    bam_files = glob.glob(args.bam)

    if not json_files:
        print("No JSON files found!")
    if not fastq_files:
        print("No FASTQ files found!")
    if not bam_files:
        print("No BAM files found!")

    fastp_data = extract_fastp_data(json_files)
    fastq_data = extract_fastq_data(fastq_files) if fastq_files else {}
    bam_data = extract_bam_data(bam_files) if bam_files else {}

    # Combine data
    all_samples = set(fastp_data.keys()) | set(fastq_data.keys()) | set(bam_data.keys())
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
            "reads_genomic": bam_data.get(sample, {}).get("reads_genomic", "NA"),
            "bases_genomic": bam_data.get(sample, {}).get("bases_genomic", "NA")
        })

    df = pd.DataFrame(summary_list)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Summary saved to {args.output}")

if __name__ == "__main__":
    main()
