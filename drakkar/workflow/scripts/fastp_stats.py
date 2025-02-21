import os
import json
import pandas as pd
import argparse
import glob

def extract_fastp_data(json_files):
    data_list = []

    for json_file in json_files:
        try:
            with open(json_file, "r") as f:
                report = json.load(f)

            sample_name = os.path.basename(json_file).replace(".json", "")
            reads_raw = report["summary"]["before_filtering"]["total_reads"]
            reads_postfiltering = report["summary"]["after_filtering"]["total_reads"]
            bases_raw = report["summary"]["before_filtering"]["total_bases"]
            bases_postfiltering = report["summary"]["after_filtering"]["total_bases"]

            data_list.append({
                "sample": sample_name,
                "reads_raw": reads_raw,
                "reads_discarded": reads_raw - reads_postfiltering,
                "bases_discarded": bases_raw - bases_postfiltering
            })

        except (KeyError, json.JSONDecodeError) as e:
            print(f"Error processing {json_file}: {e}")

    return pd.DataFrame(data_list)

def main():
    parser = argparse.ArgumentParser(description="Extract read counts from Fastp JSON reports")
    parser.add_argument("-f", "--fastp", required=True, help="Glob pattern for input JSON files (e.g., 'fastp_results/*.json')")
    parser.add_argument("-m", "--metagenome", required=False, help="Glob pattern for input JSON files (e.g., 'metagenome/*_1.fq.gz')")
    parser.add_argument("-g", "--genome", required=False, help="Glob pattern for input JSON files (e.g., 'genome/*.bam')")
    parser.add_argument("-o", "--output", required=True, help="Output filename (e.g., 'fastp_summary.tsv')")

    args = parser.parse_args()

    json_files = glob.glob(args.input)

    if not json_files:
        print("No JSON files found with the given pattern!")
        return

    df = extract_fastp_data(json_files)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Summary saved to {args.output}")

if __name__ == "__main__":
    main()
