import sys
import pandas as pd

def main():
    # Input and output file paths
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read the input file as a space-separated DataFrame
    df = pd.read_csv(
        input_file,
        sep=r"\s+",
        header=None,
        usecols=[0, 2, 4],
        names=["ko_id", "query", "evalue"],
        comment="#"
    )

    # Convert evalue to numeric for accurate comparisons (handles scientific notation)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")

    # Find the row with the lowest e-value for each query
    df_lowest_evalue = df.loc[df.groupby("query")["evalue"].idxmin()]

    # Remove annotations with too high e-values
    df_lowest_evalue = df_lowest_evalue[df_lowest_evalue["evalue"] <= 1e-5]

    # Rename columns to the desired output format
    df_lowest_evalue = df_lowest_evalue.rename(columns={"query": "gene", "ko_id": "ko_id", "evalue": "ko_e"})
    final_df = df_lowest_evalue[["gene", "ko_id", "ko_e"]]

    # Save the result to the output file
    final_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
