import argparse
import os

import pandas as pd


OUTPUT_COLUMNS = ["genome", "completeness", "contamination", "score", "size", "N50", "contig_count"]

def process_tsv_files(tsv_files, output_csv):
    combined_data = []

    for tsv_file in tsv_files:
        assembly_id = os.path.splitext(os.path.basename(tsv_file))[0]

        try:
            df = pd.read_csv(tsv_file, sep="\t")

            required_columns = {"bin_id", "completeness", "contamination"}
            if not required_columns.issubset(df.columns):
                raise ValueError(f"Missing required columns in {tsv_file}")

            df["genome"] = df["bin_id"].astype(str).apply(lambda x: f"{assembly_id}_bin_{x}.fa")

            for column in OUTPUT_COLUMNS:
                if column not in df.columns:
                    df[column] = pd.NA

            df_selected = df[OUTPUT_COLUMNS]

            combined_data.append(df_selected)

        except Exception as e:
            print(f"❌ Error processing {tsv_file}: {e}")

    if combined_data:
        final_df = pd.concat(combined_data, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=OUTPUT_COLUMNS)

    final_df.to_csv(output_csv, index=False)
    print(f"✅ Combined CSV saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and combine genome information from multiple TSV files.")
    parser.add_argument("tsv_files", nargs="+", help="List of input TSV files")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")

    args = parser.parse_args()
    process_tsv_files(args.tsv_files, args.output)
