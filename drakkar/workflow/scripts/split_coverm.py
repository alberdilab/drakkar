import argparse
import pandas as pd

def split_tsv(input_file, read_count_output, covered_bases_output):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Extract the Genome column
    genome_col = df.iloc[:, 0]

    # Create Read Count and Covered Bases dataframes
    read_counts = {'Genome': genome_col}
    covered_bases = {'Genome': genome_col}

    for col in df.columns[1:]:
        if "Read Count" in col:
            new_col_name = col.replace(" Read Count", "")
            read_counts[new_col_name] = df[col]
        elif "Covered Bases" in col:
            new_col_name = col.replace(" Covered Bases", "")
            covered_bases[new_col_name] = df[col]

    # Convert dictionaries back to DataFrame
    read_counts_df = pd.DataFrame(read_counts)
    covered_bases_df = pd.DataFrame(covered_bases)

    # Save to output files
    read_counts_df.to_csv(read_count_output, sep='\t', index=False)
    covered_bases_df.to_csv(covered_bases_output, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Split TSV file into Read Counts and Covered Bases.")
    parser.add_argument("input_file", help="Path to input TSV file")
    parser.add_argument("read_count_output", help="Path to output TSV file for Read Counts")
    parser.add_argument("covered_bases_output", help="Path to output TSV file for Covered Bases")

    args = parser.parse_args()
    split_tsv(args.input_file, args.read_count_output, args.covered_bases_output)

if __name__ == "__main__":
    main()
