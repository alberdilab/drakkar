import argparse
import pandas as pd

def create_mag_table(names, paths, output_tsv):
    """
    Creates a two-column TSV file with paths in the first column and names in the second.

    :param names: List of MAG names.
    :param paths: List of corresponding file paths.
    :param output_tsv: Path to output TSV file.
    """
    # Ensure both lists are the same length
    if len(names) != len(paths):
        raise ValueError("❌ Error: The number of names and paths must be the same.")

    # Create a DataFrame
    df = pd.DataFrame({"path": paths, "name": names})

    # Save as TSV file
    df.to_csv(output_tsv, sep="\t", index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a two-column TSV file with paths and names.")
    parser.add_argument("--names", nargs="+", required=True, help="List of MAG names.")
    parser.add_argument("--paths", nargs="+", required=True, help="List of corresponding file paths.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file.")

    args = parser.parse_args()

    # Run function
    create_mag_table(args.names, args.paths, args.output)
