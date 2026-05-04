#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import pandas as pd


def normalize_genome_name(value):
    name = Path(str(value)).name
    if name.endswith(".gz"):
        name = name[:-3]
    return name


def load_metadata(path):
    if not path:
        return None
    metadata_path = Path(path)
    if not metadata_path.exists():
        return None
    df = pd.read_csv(metadata_path)
    if "genome" not in df.columns:
        raise ValueError(f"Metadata file missing required 'genome' column: {metadata_path}")
    return df


def mean_or_na(series):
    if series is None:
        return "NA"
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return "NA"
    return f"{values.mean():.2f}"


def main():
    parser = argparse.ArgumentParser(description="Summarize input and output MAG statistics after dereplication.")
    parser.add_argument("--bins-map", required=True, help="Path to bins_to_files.json")
    parser.add_argument("--wdb", required=True, help="Path to dRep Wdb.csv")
    parser.add_argument("--ani", required=True, type=float, help="Dereplication ANI threshold")
    parser.add_argument("--metadata", required=False, help="Path to genome quality metadata CSV")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    with open(args.bins_map, "r", encoding="utf-8") as handle:
        bins_map = json.load(handle)

    input_genomes = [normalize_genome_name(path) for path in bins_map.values()]
    wdb = pd.read_csv(args.wdb)
    if "genome" not in wdb.columns:
        raise ValueError(f"Wdb missing required 'genome' column: {args.wdb}")
    output_genomes = [normalize_genome_name(value) for value in wdb["genome"].dropna().astype(str).tolist()]

    metadata = load_metadata(args.metadata)

    input_comp = input_cont = output_comp = output_cont = "NA"
    if metadata is not None:
        genome_series = metadata["genome"].astype(str).map(normalize_genome_name)
        input_df = metadata.loc[genome_series.isin(input_genomes)]
        output_df = metadata.loc[genome_series.isin(output_genomes)]
        if "completeness" in metadata.columns:
            input_comp = mean_or_na(input_df["completeness"])
            output_comp = mean_or_na(output_df["completeness"])
        if "contamination" in metadata.columns:
            input_cont = mean_or_na(input_df["contamination"])
            output_cont = mean_or_na(output_df["contamination"])

    out_df = pd.DataFrame(
        [
            {
                "input_bin_number": len(input_genomes),
                "input_bin_completeness": input_comp,
                "input_bin_contamination": input_cont,
                "dereplication_ani": f"{args.ani:.2f}",
                "output_bin_number": len(output_genomes),
                "output_bin_completeness": output_comp,
                "output_bin_contamination": output_cont,
            }
        ]
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
