import os
import pandas as pd
import argparse

def process_tsv_files(tsv_files, output_csv):
    """
    Processes multiple TSV files, extracts required columns, and saves as a CSV.

    :param tsv_files: List of input TSV file paths.
    :param output_csv: Path for the output CSV file.
    """
    combined_data = []

    for tsv_file in tsv_files:
        # Extract assembly_id from filename (e.g., SAD20.tsv -> SAD20)
        assembly_id = os.path.splitext(os.path.basename(tsv_file))[0]

        try:
            # Read TSV file
            df = pd.read_csv(tsv_file, sep="\t")

            # Check if required columns exist
            required_columns = {"bin_id", "completeness", "contamination"}
            if not required_columns.issubset(df.columns):
                raise ValueError(f"Missing required columns in {tsv_file}")

            # Create the 'genome' column
            df["genome"] = df["bin_id"].astype(str).apply(lambda x: f"{assembly_id}_bin_{x}.fa")

            # Select only the required columns
            df_selected = df[["genome", "completeness", "contamination","score","size","N50","contig_count"]]

            # Append to list
            combined_data.append(df_selected)

        except Exception as e:
            print(f"❌ Error processing {tsv_file}: {e}")

    # Concatenate all dataframes
    if combined_data:
        final_df = pd.concat(combined_data, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        print(f"✅ Combined CSV saved to {output_csv}")
    else:
        print("⚠️ No valid data found. CSV not created.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and combine genome information from multiple TSV files.")
    parser.add_argument("tsv_files", nargs="+", help="List of input TSV files")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")

    args = parser.parse_args()
    process_tsv_files(args.tsv_files, args.output)
