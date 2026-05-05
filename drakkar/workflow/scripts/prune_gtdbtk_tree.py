import argparse
import csv
from pathlib import Path

from Bio import Phylo


def load_query_genomes(batchfile_path):
    query_genomes = set()
    with open(batchfile_path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            genome_name = row[1].strip()
            if genome_name:
                query_genomes.add(genome_name)
    return query_genomes


def prune_tree(input_tree_path, batchfile_path, output_tree_path):
    query_genomes = load_query_genomes(batchfile_path)
    if not query_genomes:
        raise ValueError(f"No MAG names were found in the GTDB-Tk batch file: {batchfile_path}")

    tree = Phylo.read(input_tree_path, "newick")
    terminal_names = {terminal.name for terminal in tree.get_terminals() if terminal.name}
    matched_names = terminal_names.intersection(query_genomes)
    if not matched_names:
        raise ValueError(
            "None of the query genomes from the GTDB-Tk batch file were found in the tree tips."
        )

    for terminal_name in sorted(terminal_names - query_genomes):
        matches = list(tree.find_clades(name=terminal_name))
        if matches:
            tree.prune(matches[0])

    output_path = Path(output_tree_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output_path), "newick")


def main():
    parser = argparse.ArgumentParser(
        description="Prune GTDB-Tk classify trees to keep only the input MAG genomes."
    )
    parser.add_argument("--input-tree", required=True, help="GTDB-Tk classify tree in Newick format")
    parser.add_argument("--batchfile", required=True, help="GTDB-Tk batchfile generated from MAG names and paths")
    parser.add_argument("--output-tree", required=True, help="Output pruned Newick tree")
    args = parser.parse_args()

    prune_tree(args.input_tree, args.batchfile, args.output_tree)


if __name__ == "__main__":
    main()
