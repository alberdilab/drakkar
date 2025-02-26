import sys
import os
import pandas as pd

def parse_fasta(file_path):
    with open(file_path, "r") as fasta_file:
        contig_names = []
        for line in fasta_file:
            if line.startswith(">"):
                contig_name = line[1:].strip().split(" ")[0]  # Extracting the contig name
                contig_names.append(contig_name)
        return contig_names

def main(input_files, output_file):
    rows = []

    for file_path in input_files:
        # Infer genome name from the file name
        genome_name = os.path.basename(file_path).replace(".fna", "")

        # Extract contig names
        contig_names = parse_fasta(file_path)

        for contig_name in contig_names:
            rows.append({"genome": genome_name, "gene": contig_name})

    # Create a DataFrame and save to CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # The input files are passed as space-separated arguments
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    main(input_files, output_file)
