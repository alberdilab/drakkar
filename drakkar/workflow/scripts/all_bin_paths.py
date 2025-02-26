import os
import pandas as pd
import argparse

def extract_fasta_paths(tsv_files, output_file):
    fasta_paths = []

    for tsv_file in tsv_files:
        # Extract assembly_id from the filename (e.g., SAD20.tsv -> SAD20)
        assembly_id = os.path.splitext(os.path.basename(tsv_file))[0]

        # Read the TSV file and extract bin IDs
        try:
            df = pd.read_csv(tsv_file, sep="\t")
            if "bin_id" not in df.columns:
                raise ValueError(f"Column 'bin_id' not found in {tsv_file}")

            bin_ids = df["bin_id"].unique()

            # Generate paths
            for bin_id in bin_ids:
                fasta_path = f"cataloging/binette/{assembly_id}/final_bins/{assembly_id}_bin_{bin_id}.fa"
                fasta_paths.append(fasta_path)

        except Exception as e:
            print(f"❌ Error processing {tsv_file}: {e}")

    # Write output file
    with open(output_file, "w") as f:
        for path in fasta_paths:
            f.write(path + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate FASTA file paths from TSV bin_id columns.")
    parser.add_argument("tsv_files", nargs="+", help="List of TSV files to process")
    parser.add_argument("-o", "--output", required=True, help="Output text file to store paths")

    args = parser.parse_args()
    extract_fasta_paths(args.tsv_files, args.output)
