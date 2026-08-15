"""Helpers to read Binette final bin quality reports.

Binette >= 1.2 renamed the report column `bin_id` to `name` and writes the bin
FASTA files as `final_bins/<name>.fa`, where the name is built as
`<prefix>_bin<n>`. DRAKKAR keeps plain numeric bin ids in its own file names
(`<assembly>_bin_<n>.fa`), so the numeric part is parsed back out of the name.
Reports written by Binette <= 1.1 (`bin_id` column) are still accepted.
"""

import re

import pandas as pd

# Prefix passed to binette --prefix in the cataloging workflow.
BINETTE_BIN_PREFIX = "binette"

# Columns written by binette >= 1.2 in final_bins_quality_reports.tsv.
REPORT_COLUMNS = [
    "name",
    "origin",
    "is_original",
    "original_name",
    "completeness",
    "contamination",
    "score",
    "checkm2_model",
    "size",
    "N50",
    "coding_density",
    "contig_count",
]

_BIN_NAME_PATTERN = re.compile(r"bin(\d+)$")


def bin_id_from_name(name):
    """Return the numeric bin id encoded in a binette bin name."""
    match = _BIN_NAME_PATTERN.search(str(name))
    return match.group(1) if match else str(name)


def bin_id_column(table):
    """Return the bin id of every row of a binette final bin quality report."""
    if "bin_id" in table.columns:
        values = table["bin_id"]
    elif "name" in table.columns:
        values = table["name"]
    else:
        return []
    return ["" if pd.isna(value) else bin_id_from_name(value) for value in values]


def bin_ids_from_table(table):
    """Return the unique bin ids of a binette final bin quality report."""
    return list(dict.fromkeys(bin_id for bin_id in bin_id_column(table) if bin_id))
