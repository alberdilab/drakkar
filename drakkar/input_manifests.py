import json
import os
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from drakkar.downloads import (
    _normalized_value,
    _validate_paired_read_maps,
    download_to_cache,
    is_url,
    resolve_input_manifest,
    resolve_preprocessed_read_lists,
    resolve_sample_read_lists,
)
from drakkar.input_errors import InputFileError, report_input_resolution_errors, require_non_empty_file
from drakkar.output import print

ASSEMBLY_COLUMN_CANDIDATES = ("assembly", "coassembly")

def check_reference_columns(file_path):
    """Checks if a file contains 'reference_name' and 'reference_path' columns with non-NA values."""
    # Read the file (assumed to be TSV, change sep="," for CSV)
    df = pd.read_csv(file_path, sep="\t")
    # Check if required columns exist
    required_columns = {"reference_name", "reference_path"}
    if not required_columns.issubset(df.columns):
        return False
    # Check if both columns have at least one non-NA value
    if df["reference_name"].dropna().empty or df["reference_path"].dropna().empty:
        return False
    return True

def get_assembly_column_name(df):
    for column in ASSEMBLY_COLUMN_CANDIDATES:
        if column in df.columns:
            return column
    return None

def check_assembly_column(file_path):
    """Checks if a file contains the preferred 'assembly' column or legacy 'coassembly' values."""
    # Read the file (assumed to be TSV, change sep="," for CSV)
    df = pd.read_csv(file_path, sep="\t")
    assembly_column = get_assembly_column_name(df)
    if not assembly_column:
        return False
    values = df[assembly_column].dropna().astype(str).str.strip()
    if values.empty or (values == "").all():
        return False
    return True

def file_samples_to_json(infofile, output):
    df = pd.read_csv(infofile, sep="\t")

    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)
    errors = []

    for idx, row in df.iterrows():
        try:
            sample_name, read1_paths, read2_paths = resolve_sample_read_lists(row, idx + 1, output)
            SAMPLE_TO_READS1[sample_name].extend(read1_paths)
            SAMPLE_TO_READS2[sample_name].extend(read2_paths)
        except InputFileError as exc:
            errors.append(str(exc))

    if errors:
        report_input_resolution_errors(errors)
    if not SAMPLE_TO_READS1 and not SAMPLE_TO_READS2:
        report_input_resolution_errors([f"No sample rows were found in {infofile}."])

    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def file_references_to_json(infofile, output):
    df = pd.read_csv(infofile, sep="\t")

    required_columns = {"reference_name", "reference_path"}
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        report_input_resolution_errors([
            f"Reference table {infofile} is missing required column(s): {', '.join(missing_columns)}"
        ])

    REFERENCE_TO_FILE = {}
    errors = []
    for idx, row in df.iterrows():
        row_number = idx + 1
        ref_name_value = _normalized_value(row.get("reference_name"))
        ref_path_value = _normalized_value(row.get("reference_path"))
        if not ref_name_value:
            errors.append(f"Missing reference_name on row {row_number} of {infofile}.")
            continue
        if not ref_path_value:
            errors.append(f"Missing reference_path on row {row_number} of {infofile}.")
            continue
        if is_url(ref_path_value):
            try:
                resolved_ref_path = download_to_cache(
                    ref_path_value,
                    ref_name_value,
                    "reference_path",
                    output,
                    cache_subdir="references_cache",
                )
            except InputFileError as exc:
                errors.append(str(exc))
                continue
        else:
            resolved_ref_path = str(Path(ref_path_value).resolve())
            try:
                require_non_empty_file(resolved_ref_path, f"reference_path on row {row_number}")
            except InputFileError as exc:
                errors.append(str(exc))
                continue
        REFERENCE_TO_FILE[ref_name_value] = resolved_ref_path

    if errors:
        report_input_resolution_errors(errors)
    if not REFERENCE_TO_FILE:
        report_input_resolution_errors([f"No reference entries were found in {infofile}."])

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    SAMPLE_TO_REFERENCE = dict(zip(df["sample"], df["reference_name"]))
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

def argument_samples_to_json(argument, output):
    # Define the directory containing the raw reads
    READS_DIR = Path(argument).resolve()

    # Initialize dictionaries
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(READS_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(READS_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    SAMPLE_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    SAMPLE_TO_READS2[sample_name].append(full_path)

    _validate_paired_read_maps(SAMPLE_TO_READS1, SAMPLE_TO_READS2, READS_DIR, "input directory")

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_references_to_json(argument, sample_to_reads, output):
    try:
        if is_url(argument):
            reference_path = download_to_cache(
                argument,
                "reference",
                "reference",
                output,
                cache_subdir="references_cache",
            )
        else:
            reference_path = str(Path(argument).resolve())
            require_non_empty_file(reference_path, "reference argument")
    except InputFileError as exc:
        report_input_resolution_errors([str(exc)])

    REFERENCE_TO_FILE = {"reference": [reference_path]}
    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    with open(f"{sample_to_reads}", "r") as f:
        SAMPLE_TO_READS = json.load(f)

    SAMPLE_TO_REFERENCE = {sample: "reference" for sample in SAMPLE_TO_READS.keys()}
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

def file_preprocessed_to_json(infofile, output):
    df = pd.read_csv(infofile, sep="\t")

    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)
    errors = []

    for idx, row in df.iterrows():
        try:
            sample_name, read1_paths, read2_paths = resolve_preprocessed_read_lists(row, idx + 1, output)
            SAMPLE_TO_READS1[sample_name].extend(read1_paths)
            SAMPLE_TO_READS2[sample_name].extend(read2_paths)
        except InputFileError as exc:
            errors.append(str(exc))

    if errors:
        report_input_resolution_errors(errors)
    if not SAMPLE_TO_READS1 and not SAMPLE_TO_READS2:
        report_input_resolution_errors([f"No sample rows were found in {infofile}."])

    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/preprocessed_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/preprocessed_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_preprocessed_to_json(argument, output):
    # Define the directory containing the raw reads
    PREPROCESSED_DIR = Path(argument).resolve()

    # Initialize dictionaries
    PREPROCESSED_TO_READS1 = defaultdict(list)
    PREPROCESSED_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(PREPROCESSED_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(PREPROCESSED_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    PREPROCESSED_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    PREPROCESSED_TO_READS2[sample_name].append(full_path)

    _validate_paired_read_maps(PREPROCESSED_TO_READS1, PREPROCESSED_TO_READS2, PREPROCESSED_DIR, "preprocessed input directory")

    # Convert defaultdict to standard dict (optional)
    PREPROCESSED_TO_READS1 = dict(PREPROCESSED_TO_READS1)
    PREPROCESSED_TO_READS2 = dict(PREPROCESSED_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/preprocessed_to_reads1.json", "w") as f:
        json.dump(PREPROCESSED_TO_READS1, f)

    with open(f"{output}/data/preprocessed_to_reads2.json", "w") as f:
        json.dump(PREPROCESSED_TO_READS2, f)

def file_assemblies_to_json(infofile=None, samples=None, individual=False, all=False, output=False):

    assemblies = defaultdict(list)

    if infofile is not None:
        df = pd.read_csv(infofile, sep="\t")
        assembly_column = get_assembly_column_name(df)
        if assembly_column:
            for _, row in df.iterrows():
                sample = row["sample"]
                assembly_value = row[assembly_column]
                if pd.isna(assembly_value) or str(assembly_value).strip() == "":
                    continue
                assembly_list = str(assembly_value).split(",")

                for assembly in assembly_list:
                    assembly = assembly.strip()
                    if assembly:
                        assemblies[assembly].append(sample)

    if samples:
        if individual:
            for sample in samples:
                assemblies[sample].append(sample)

        if all:
            assemblies["all"].extend(samples)

    ASSEMBLY_TO_SAMPLE = {key: sorted(set(value)) for key, value in assemblies.items()}

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/assembly_to_samples.json", "w") as f:
        json.dump(ASSEMBLY_TO_SAMPLE, f)

def file_coverages_to_json(infofile=None, samples=None, output=False):

    coverage_groups = defaultdict(list)
    has_coverage_column = False

    if infofile is not None:
        df = pd.read_csv(infofile, sep="\t")
        if "coverage" in df.columns:
            has_coverage_column = True
            for _, row in df.iterrows():
                sample = row["sample"]
                coverage_value = row["coverage"]
                if pd.isna(coverage_value) or str(coverage_value).strip() == "":
                    coverage_value = sample
                for group in str(coverage_value).split(","):
                    group = group.strip()
                    if group:
                        coverage_groups[group].append(sample)

    if not has_coverage_column and samples:
        coverage_groups["all"].extend(samples)

    COVERAGE_TO_SAMPLES = {
        key: sorted(set(value)) for key, value in coverage_groups.items()
    }

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/coverage_to_samples.json", "w") as f:
        json.dump(COVERAGE_TO_SAMPLES, f)
    return has_coverage_column

def file_bins_to_json(paths_file=None, output=False):
    fasta_dict = {}

    try:
        paths_file = resolve_input_manifest(paths_file, output, label="bins_file")
        require_non_empty_file(paths_file, "bins_file manifest")
    except InputFileError as exc:
        report_input_resolution_errors([str(exc)])

    if not os.path.isfile(paths_file):
        report_input_resolution_errors([f"Bin file not found: {paths_file}"])

    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)
    errors = []

    with open(paths_file, "r") as f:
        for line_number, line in enumerate(f, start=1):
            full_path = line.strip()
            if not full_path:
                continue
            if is_url(full_path):
                try:
                    full_path = download_to_cache(
                        full_path,
                        "",
                        "genome",
                        output,
                        cache_subdir="genomes_cache",
                        preserve_basename=True,
                    )
                except InputFileError as exc:
                    errors.append(str(exc))
                    continue
            else:
                try:
                    require_non_empty_file(full_path, f"genome file listed in {paths_file}:{line_number}")
                except InputFileError as exc:
                    errors.append(str(exc))
                    continue

            filename = fasta_re.sub("", os.path.basename(full_path))
            fasta_dict[filename] = full_path

    if not fasta_dict and not errors:
        errors.append(f"No genome FASTA paths were found in {paths_file}.")

    if errors:
        report_input_resolution_errors(errors)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def path_bins_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        report_input_resolution_errors([f"Folder not found: {folder_path}"])

    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)
    errors = []

    # Iterate over all files in the folder
    for file_name in os.listdir(folder_path):
        if fasta_re.search(file_name):
            full_path = os.path.join(folder_path, file_name)
            try:
                require_non_empty_file(full_path, f"bin FASTA file in {folder_path}")
            except InputFileError as exc:
                errors.append(str(exc))
                continue
            file_id = fasta_re.sub("", file_name)
            fasta_dict[file_id] = full_path

    if not fasta_dict and not errors:
        errors.append(f"No FASTA files were found in {folder_path}.")
    if errors:
        report_input_resolution_errors(errors)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def file_mags_to_json(paths_file=None, output=False):
    fasta_dict = {}
    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)

    try:
        paths_file = resolve_input_manifest(paths_file, output, label="mags_file")
        require_non_empty_file(paths_file, "mags_file manifest")
    except InputFileError as exc:
        report_input_resolution_errors([str(exc)])

    if not os.path.isfile(paths_file):
        report_input_resolution_errors([f"MAG file not found: {paths_file}"])

    errors = []

    with open(paths_file, "r") as f:
        for line_number, line in enumerate(f, start=1):
            full_path = line.strip()
            if not full_path:
                continue
            if is_url(full_path):
                try:
                    full_path = download_to_cache(
                        full_path,
                        "",
                        "genome",
                        output,
                        cache_subdir="genomes_cache",
                        preserve_basename=True,
                    )
                except InputFileError as exc:
                    errors.append(str(exc))
                    continue
            else:
                try:
                    require_non_empty_file(full_path, f"genome file listed in {paths_file}:{line_number}")
                except InputFileError as exc:
                    errors.append(str(exc))
                    continue

            filename = fasta_re.sub("", os.path.basename(full_path))
            fasta_dict[filename] = full_path

    if not fasta_dict and not errors:
        errors.append(f"No genome FASTA paths were found in {paths_file}.")

    if errors:
        report_input_resolution_errors(errors)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/mags_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def path_mags_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        report_input_resolution_errors([f"Folder not found: {folder_path}"])

    # Compile a regex that matches .fa/.fna/.fasta, optionally followed by .gz
    FASTA_RE = re.compile(r'\.(?:fa|fna|fasta)(?:\.gz)?$', re.IGNORECASE)

    # Iterate over all files in the folder
    fasta_dict = {}
    errors = []
    for fname in os.listdir(folder_path):
        if FASTA_RE.search(fname):
            full_path = os.path.join(folder_path, fname)
            try:
                require_non_empty_file(full_path, f"MAG FASTA file in {folder_path}")
            except InputFileError as exc:
                errors.append(str(exc))
                continue
            sample_id = FASTA_RE.sub('', fname)
            fasta_dict[sample_id] = full_path

    if not fasta_dict and not errors:
        errors.append(f"No FASTA files were found in {folder_path}.")
    if errors:
        report_input_resolution_errors(errors)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/mags_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def microdiversity_selection_to_json(coverage, mincov=0.5, minsamp=10):
    df = pd.read_csv(csv_input, sep='\t', index_col=0)
    genome_sample_dict = {}
    for genome, row in df.iterrows():
        selected = row[row > mincov].index.tolist()
        if len(selected) >= minsamp:
            genome_sample_dict[genome] = selected

    with open(f"{output}/data/microdiversity_selection.json", "w") as f:
        json.dump(genome_sample_dict, f, indent=4)

def preprocessing_summary(summary_table, bar_width=50):
    """
    Prints a single horizontal stacked barplot representing the AVERAGE percentage of
    bases_discarded, bases_host, and bases_metagenomic across all samples.

    Uses:
    - '░' (light block) for discarded bases
    - '▒' (medium block) for host bases
    - '█' (full block) for metagenomic bases
    - '│' (vertical separator) between categories
    """

    df = pd.read_csv(summary_table, sep="\t")

    # Compute total sums
    total_discarded = df["bases_discarded"].sum()
    total_host = df["bases_host"].sum()
    total_metagenomic = df["bases_metagenomic"].sum()

    total_bases = total_discarded + total_host + total_metagenomic
    if total_bases == 0:
        print("No data available")
        return

    # Compute AVERAGE percentages
    pct_discarded = (total_discarded / total_bases) * 100
    pct_host = (total_host / total_bases) * 100
    pct_metagenomic = (total_metagenomic / total_bases) * 100

    # Compute character count for each section
    discarded_chars = round((pct_discarded / 100) * bar_width)
    host_chars = round((pct_host / 100) * bar_width)
    metagenomic_chars = bar_width - discarded_chars - host_chars  # Ensure total width matches

    # Construct the visual bar
    bar = "░" * discarded_chars + "│" + "▒" * host_chars + "│" + "█" * metagenomic_chars

    # Print the averages and the stacked barplot
    print(f"░ Discarded: {pct_discarded:.1f}%, ▒ Host: {pct_host:.1f}%, █ Metagenomic: {pct_metagenomic:.1f}%")
    print(f"│{bar}│")

def file_transcriptome_to_json(infofile, output):
    df = pd.read_csv(infofile, sep="\t")

    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)
    errors = []

    for idx, row in df.iterrows():
        try:
            sample_name, read1_paths, read2_paths = resolve_sample_read_lists(row, idx + 1, output)
            SAMPLE_TO_READS1[sample_name].extend(read1_paths)
            SAMPLE_TO_READS2[sample_name].extend(read2_paths)
        except InputFileError as exc:
            errors.append(str(exc))

    if errors:
        report_input_resolution_errors(errors)
    if not SAMPLE_TO_READS1 and not SAMPLE_TO_READS2:
        report_input_resolution_errors([f"No sample rows were found in {infofile}."])

    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/transcriptome_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/transcriptome_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_transcriptome_to_json(argument, output):
    # Define the directory containing the raw reads
    TRANSCRIPTOME_DIR = Path(argument).resolve()

    # Initialize dictionaries
    TRANSCRIPTOME_TO_READS1 = defaultdict(list)
    TRANSCRIPTOME_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(TRANSCRIPTOME_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(TRANSCRIPTOME_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    TRANSCRIPTOME_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    TRANSCRIPTOME_TO_READS2[sample_name].append(full_path)

    _validate_paired_read_maps(TRANSCRIPTOME_TO_READS1, TRANSCRIPTOME_TO_READS2, TRANSCRIPTOME_DIR, "transcriptome input directory")

    # Convert defaultdict to standard dict (optional)
    TRANSCRIPTOME_TO_READS1 = dict(TRANSCRIPTOME_TO_READS1)
    TRANSCRIPTOME_TO_READS2 = dict(TRANSCRIPTOME_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/transcriptome_to_reads1.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS1, f)

    with open(f"{output}/data/transcriptome_to_reads2.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS2, f)
