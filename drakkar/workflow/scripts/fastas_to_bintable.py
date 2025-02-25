import os
import argparse
import pandas as pd

def parse_fasta_files(directory, extension, output_file):
    data = []

    # Iterate through files in directory
    for filename in os.listdir(directory):
        if filename.endswith(extension):
            file_path = os.path.join(directory, filename)
            file_name_without_ext = os.path.splitext(filename)[0]

            # Read the fasta file manually
            with open(file_path, "r") as fasta_file:
                for line in fasta_file:
                    if line.startswith(">"):
                        contig_name = line.strip().lstrip(">")
                        data.append([contig_name, file_name_without_ext])

    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(data, columns=["Contig_Name", "File_Name"])
    df.to_csv(output_file, index=False, header=False, sep="\t")
    print(f"Table saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract contig names from FASTA files and save as a table.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing FASTA files")
    parser.add_argument("-e", "--extension", required=True, help="FASTA file extension (e.g., .fa, .fasta, .fna)")
    parser.add_argument("-o", "--output", required=True, help="Output tsv file")

    args = parser.parse_args()

    parse_fasta_files(args.directory, args.extension, args.output)
