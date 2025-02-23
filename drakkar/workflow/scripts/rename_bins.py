import os
import sys

def rename_fasta_headers(assembly_name, input_file, output_file):
    # Extract original file name and generate new file name
    input_basename = os.path.basename(input_file)
    new_filename = f"{assembly_name}_{input_basename}"  # Project-prefixed filename

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                # Modify the header: >project@original_header
                new_header = f">{assembly_name}@{line[1:].strip()}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)  # Keep sequence lines unchanged

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_fasta.py <assembly_name> <input_fasta> <output_fasta>")
        sys.exit(1)

    # Read command-line arguments
    assembly_name = sys.argv[1]
    input_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    # Run renaming function
    rename_fasta_headers(assembly_name, input_fasta, output_fasta)
