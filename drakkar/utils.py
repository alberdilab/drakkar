import csv
import os
import sys
import yaml
import re
import json
import time
import pandas as pd
import shutil
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen
from pathlib import Path
from collections import defaultdict

from drakkar import __version__
from drakkar.ascii import (
    DRAKKAR_LOGO_ART,
    DRAKKAR_SHIP_ART,
    END_ART,
    INTRO_ROWS,
    UNLOCK_ART,
    UPDATE_SUCCESS_ASCII_TEMPLATE,
)
from drakkar.output import Text as RichText, print, prompt

class DownloadError(Exception):
    pass

DRAKKAR_SHIP_STYLE = "bold #5f9ea0"
DRAKKAR_LOGO_STYLE = "bold #d6a642"
DRAKKAR_INTRO_STYLE = "bold #b7c7d3"
DRAKKAR_VERSION_BADGE_STYLE = DRAKKAR_INTRO_STYLE
BANNER_DELAY_SECONDS = 0.3
ASSEMBLY_COLUMN_CANDIDATES = ("assembly", "coassembly")
READ1_BASENAME_PATTERN = re.compile(r"(?:^|[._-])(?:R?1)(?:[._-]|$)", re.IGNORECASE)
READ2_BASENAME_PATTERN = re.compile(r"(?:^|[._-])(?:R?2)(?:[._-]|$)", re.IGNORECASE)

def is_url(value):
    parsed = urlparse(str(value))
    return parsed.scheme in {"http", "https", "ftp"}

def _has_value(value):
    return not (pd.isna(value) or str(value).strip() == "")

def _normalized_value(value):
    if pd.isna(value):
        return ""
    return str(value).strip()

def download_to_cache(url, sample_name, column_name, output, cache_subdir="reads_cache", preserve_basename=False, max_retries=3):
    cache_dir = os.path.join(output, "data", cache_subdir)
    os.makedirs(cache_dir, exist_ok=True)
    parsed = urlparse(url)
    basename = os.path.basename(parsed.path) or f"{sample_name}_{column_name}"
    if preserve_basename or not sample_name:
        filename = basename
    else:
        filename = f"{sample_name}_{basename}"
    dest_path = os.path.join(cache_dir, filename)
    if os.path.exists(dest_path):
        print(f"Using cached {column_name} for {sample_name}: {dest_path}", flush=True)
        return dest_path

    print(f"Downloading {column_name} for {sample_name} from {url}", flush=True)
    tmp_path = f"{dest_path}.tmp"
    for attempt in range(1, max_retries + 1):
        try:
            with urlopen(url) as response, open(tmp_path, "wb") as handle:
                shutil.copyfileobj(response, handle)
            os.replace(tmp_path, dest_path)
            print(f"Saved {column_name} for {sample_name} to {dest_path}", flush=True)
            return dest_path
        except Exception as exc:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            if attempt < max_retries:
                delay = 5 * (3 ** (attempt - 1))
                print(f"WARNING: Download attempt {attempt}/{max_retries} failed for {url}: {exc}. Retrying in {delay}s...", flush=True)
                time.sleep(delay)
            else:
                raise DownloadError(f"Failed to download {url} after {max_retries} attempts: {exc}")

def resolve_input_manifest(manifest_path, output, cache_subdir="manifests_cache", label="manifest"):
    if is_url(manifest_path):
        if not output:
            raise ValueError(f"An output directory is required to download the remote {label}.")
        return download_to_cache(
            manifest_path,
            "",
            label,
            output,
            cache_subdir=cache_subdir,
            preserve_basename=True,
        )
    return manifest_path

def _resolve_reads_path(read_value, sample_name, column_name, output, row_number):
    if is_url(read_value):
        return download_to_cache(read_value, sample_name, column_name, output)

    resolved_path = str(Path(read_value).resolve())
    if not os.path.exists(resolved_path):
        print(f"ERROR: {column_name} file not found on row {row_number}: {resolved_path}")
        sys.exit(1)
    return resolved_path

def _normalize_ena_fastq_url(url):
    normalized = _normalized_value(url)
    if not normalized:
        return None
    if is_url(normalized):
        return normalized
    return f"https://{normalized.lstrip('/')}"

def _split_paired_fastq_urls(urls, accession):
    read1_urls = []
    read2_urls = []
    unclassified = []

    for url in urls:
        basename = os.path.basename(urlparse(url).path)
        if READ1_BASENAME_PATTERN.search(basename):
            read1_urls.append(url)
        elif READ2_BASENAME_PATTERN.search(basename):
            read2_urls.append(url)
        else:
            unclassified.append(url)

    if unclassified:
        if len(urls) == 2 and not read1_urls and not read2_urls:
            return [urls[0]], [urls[1]]
        raise DownloadError(
            "Could not determine forward and reverse FASTQ files "
            f"for accession {accession}. Use explicit rawreads1/rawreads2 paths instead."
        )

    if not read1_urls or not read2_urls or len(read1_urls) != len(read2_urls):
        raise DownloadError(
            f"Accession {accession} did not resolve to a balanced paired-end FASTQ set. "
            "Use explicit rawreads1/rawreads2 paths instead."
        )

    return sorted(read1_urls), sorted(read2_urls)

def resolve_accession_to_reads(accession, sample_name, output, row_number):
    accession_value = _normalized_value(accession)
    query_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        + urlencode(
            {
                "accession": accession_value,
                "result": "read_run",
                "fields": "run_accession,library_layout,fastq_ftp",
                "format": "tsv",
            }
        )
    )
    payload = None
    for attempt in range(1, 4):
        try:
            with urlopen(query_url) as response:
                payload = response.read().decode("utf-8")
            break
        except Exception as exc:
            if attempt < 3:
                delay = 5 * (3 ** (attempt - 1))
                print(f"WARNING: ENA API query attempt {attempt}/3 failed for accession {accession_value}: {exc}. Retrying in {delay}s...", flush=True)
                time.sleep(delay)
            else:
                raise DownloadError(f"Failed to query ENA API for accession {accession_value} on row {row_number}: {exc}")

    reader = csv.DictReader(payload.splitlines(), delimiter="\t")
    rows = [row for row in reader if any(_normalized_value(value) for value in row.values())]
    if not rows:
        raise DownloadError(
            f"No ENA/SRA FASTQ metadata were found for accession {accession_value} "
            f"on row {row_number}."
        )

    matching_rows = [row for row in rows if _normalized_value(row.get("run_accession")) == accession_value]
    selected_row = matching_rows[0] if matching_rows else rows[0]
    library_layout = _normalized_value(selected_row.get("library_layout")).upper()
    fastq_field = _normalized_value(selected_row.get("fastq_ftp"))

    if library_layout and library_layout != "PAIRED":
        raise DownloadError(
            f"Accession {accession_value} on row {row_number} is reported as "
            f"{library_layout.lower()}, but DRAKKAR expects paired-end reads here."
        )

    if not fastq_field:
        raise DownloadError(
            f"No FASTQ download links were found for accession {accession_value} "
            f"on row {row_number}."
        )

    fastq_urls = [_normalize_ena_fastq_url(item) for item in fastq_field.split(";")]
    fastq_urls = [url for url in fastq_urls if url]
    if not fastq_urls:
        raise DownloadError(
            f"Accession {accession_value} on row {row_number} returned empty FASTQ links."
        )

    forward_urls, reverse_urls = _split_paired_fastq_urls(fastq_urls, accession_value)
    reads1 = [
        download_to_cache(url, sample_name, "rawreads1", output)
        for url in forward_urls
    ]
    reads2 = [
        download_to_cache(url, sample_name, "rawreads2", output)
        for url in reverse_urls
    ]
    return reads1, reads2

def resolve_sample_read_lists(row, row_number, output):
    sample = row.get("sample")
    rawreads1 = row.get("rawreads1")
    rawreads2 = row.get("rawreads2")
    accession = row.get("accession")

    if not _has_value(sample):
        print(f"ERROR: Missing value in column 'sample' on row {row_number} of the sample info file.")
        sys.exit(1)

    sample_name = _normalized_value(sample)
    accession_value = _normalized_value(accession)
    rawreads1_value = _normalized_value(rawreads1)
    rawreads2_value = _normalized_value(rawreads2)

    if accession_value:
        if rawreads1_value or rawreads2_value:
            print(
                f"ERROR: Row {row_number} of the sample info file mixes 'accession' with "
                "'rawreads1'/'rawreads2'. Use either accession or explicit read files."
            )
            sys.exit(1)
        reads1, reads2 = resolve_accession_to_reads(accession_value, sample_name, output, row_number)
        return sample_name, reads1, reads2

    for column, value in (("rawreads1", rawreads1), ("rawreads2", rawreads2)):
        if not _has_value(value):
            print(f"ERROR: Missing value in column '{column}' on row {row_number} of the sample info file.")
            sys.exit(1)

    reads1 = [_resolve_reads_path(rawreads1_value, sample_name, "rawreads1", output, row_number)]
    reads2 = [_resolve_reads_path(rawreads2_value, sample_name, "rawreads2", output, row_number)]
    return sample_name, reads1, reads2

def _version_badge_lines(version):
    label = f"v{version}"
    width = len(label) + 2
    return [
        "╭" + "─" * width + "╮",
        f"│ {label} │",
        "╰" + "─" * width + "╯",
    ]

def _ascii_block_with_bottom_right_badge(ascii_text, version):
    lines = ascii_text.splitlines()
    if not lines:
        return ascii_text

    non_empty_line_indexes = [index for index, line in enumerate(lines) if line.strip()]
    if not non_empty_line_indexes:
        return ascii_text

    badge_lines = _version_badge_lines(version)
    badge_width = max(len(line) for line in badge_lines)
    block_width = max(len(line) for line in lines)
    badge_start = max(0, block_width - badge_width)
    first_badge_line = max(
        non_empty_line_indexes[0],
        non_empty_line_indexes[-1] - len(badge_lines) + 1,
    )
    target_indexes = [first_badge_line + offset for offset in range(len(badge_lines))]
    content_width = max(
        len(lines[index].rstrip()) for index in target_indexes if index < len(lines)
    )
    if content_width > badge_start:
        badge_start = content_width + 2

    for offset, badge_line in enumerate(badge_lines):
        line_index = first_badge_line + offset
        if line_index >= len(lines):
            lines.append("")
        lines[line_index] = lines[line_index].ljust(badge_start) + badge_line

    return "\n".join(lines)

def _ascii_block_width(ascii_text):
    return max((len(line) for line in ascii_text.splitlines()), default=0)

def _center_block(block, target_width):
    block_width = _ascii_block_width(block)
    padding = max(0, (target_width - block_width) // 2)
    return "\n".join(
        (" " * padding + line) if line.strip() else line
        for line in block.splitlines()
    )

def _intro_box(target_width):
    content_width = max(len(row) for row in INTRO_ROWS)
    inner_width = content_width + 4
    lines = [
        "╭" + "─" * inner_width + "╮",
        *[f"│  {row.center(content_width)}  │" for row in INTRO_ROWS],
        "╰" + "─" * inner_width + "╯",
    ]
    return _center_block("\n".join(lines), target_width)

def _styled_ascii_art(ascii_text, base_style, accent_chunks=None, accent_style=None):
    if RichText is None:
        return ascii_text

    rendered = RichText(ascii_text, style=base_style, no_wrap=True, overflow="ignore")
    if accent_chunks and accent_style:
        search_start = 0
        for chunk in accent_chunks:
            start = ascii_text.find(chunk, search_start)
            if start == -1:
                start = ascii_text.find(chunk)
            if start == -1:
                continue
            rendered.stylize(accent_style, start, start + len(chunk))
            search_start = start + len(chunk)
    return rendered

def get_drakkar_banner_blocks(include_intro=True):
    ascii_ship = _ascii_block_with_bottom_right_badge(DRAKKAR_SHIP_ART, __version__)
    blocks = [
        (ascii_ship, DRAKKAR_SHIP_STYLE),
        (DRAKKAR_LOGO_ART, DRAKKAR_LOGO_STYLE),
    ]
    if include_intro:
        ascii_intro = _intro_box(_ascii_block_width(DRAKKAR_LOGO_ART))
        blocks.append((ascii_intro, DRAKKAR_INTRO_STYLE))
    return blocks

def get_drakkar_banner_renderables(include_intro=True):
    ascii_ship = _ascii_block_with_bottom_right_badge(DRAKKAR_SHIP_ART, __version__)
    renderables = [
        _styled_ascii_art(
            ascii_ship,
            DRAKKAR_SHIP_STYLE,
            accent_chunks=_version_badge_lines(__version__),
            accent_style=DRAKKAR_VERSION_BADGE_STYLE,
        ),
        _styled_ascii_art(DRAKKAR_LOGO_ART, DRAKKAR_LOGO_STYLE),
    ]
    if include_intro:
        ascii_intro = _intro_box(_ascii_block_width(DRAKKAR_LOGO_ART))
        renderables.append(_styled_ascii_art(ascii_intro, DRAKKAR_INTRO_STYLE))
    return renderables

def _banner_animation_enabled(stream=None):
    if os.environ.get("DRAKKAR_NO_ANIMATION"):
        return False
    stream = stream or sys.stdout
    return bool(getattr(stream, "isatty", lambda: False)())

def display_banner_sequence(renderables, file=None, delay_after=False):
    delay_enabled = _banner_animation_enabled(file)
    rendered = list(renderables)
    for index, renderable in enumerate(rendered):
        is_last = index == len(rendered) - 1
        print(renderable, file=file)
        if delay_enabled and (not is_last or delay_after):
            time.sleep(BANNER_DELAY_SECONDS)

def display_drakkar(file=None):
    display_banner_sequence(get_drakkar_banner_renderables(), file=file)

def display_unlock():
    print(UNLOCK_ART)

def display_end():
    print(END_ART)


def _format_update_success_ascii(version):
    return UPDATE_SUCCESS_ASCII_TEMPLATE.replace("X.X.XX", str(version).center(len("X.X.XX")), 1)


def display_update_success(version):
    print(_styled_ascii_art(_format_update_success_ascii(version), DRAKKAR_SHIP_STYLE))

def is_snakemake_locked(workdir: str) -> bool:
    locks_dir = os.path.join(workdir, ".snakemake", "locks")
    return os.path.isdir(locks_dir) and len(os.listdir(locks_dir)) > 0

def check_screen_session():
    """Checks if the script is running inside a 'screen' session. If not, warns the user and asks for confirmation."""
    if "STY" not in os.environ:
        print("\n ⚠️   WARNING: You are not running this script inside a 'screen' session.")
        print("     Running long processes outside of 'screen' may cause issues if your session is disconnected.")
        if RichText is None:
            print("\n 📌   To start a screen session, use: `screen -S mysession`")
        else:
            screen_message = RichText("\n 📌   To start a screen session, use: ")
            screen_message.append("`screen -S mysession`", style="drakkar.code")
            print(screen_message)
        print("         Then run this script inside the screen session.\n")

        # Prompt user to continue
        while True:
            user_input = prompt("    👉 Type '1' to ignore this warning and continue, or Ctrl+C to exit: ")
            if user_input.strip() == "1":
                break  # Continue execution
            else:
                print("     Invalid input. Please type '1' to continue.")

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
        except DownloadError as exc:
            errors.append(str(exc))

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)


def file_references_to_json(infofile, output):
    df = pd.read_csv(infofile, sep="\t")

    REFERENCE_TO_FILE = {}
    errors = []
    for ref_name, ref_path in zip(df["reference_name"], df["reference_path"]):
        ref_name_value = str(ref_name)
        ref_path_value = str(ref_path)
        if is_url(ref_path_value):
            try:
                resolved_ref_path = download_to_cache(
                    ref_path_value,
                    ref_name_value,
                    "reference_path",
                    output,
                    cache_subdir="references_cache",
                )
            except DownloadError as exc:
                errors.append(str(exc))
                continue
        else:
            resolved_ref_path = str(Path(ref_path_value).resolve())
        REFERENCE_TO_FILE[ref_name_value] = resolved_ref_path

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    SAMPLE_TO_REFERENCE = dict(zip(df["sample"], df["reference_name"]))
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

#def file_assemblies_to_json(infofile):

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

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_references_to_json(argument, sample_to_reads, output):
    if is_url(argument):
        reference_path = download_to_cache(
            argument,
            "reference",
            "reference",
            output,
            cache_subdir="references_cache",
        )
    else:
        reference_path = argument

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
            sample_name, read1_paths, read2_paths = resolve_sample_read_lists(row, idx + 1, output)
            SAMPLE_TO_READS1[sample_name].extend(read1_paths)
            SAMPLE_TO_READS2[sample_name].extend(read2_paths)
        except DownloadError as exc:
            errors.append(str(exc))

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

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

# Create dictionary of bin names and paths from the input file
def file_bins_to_json(paths_file=None, output=False):
    fasta_dict = {}

    paths_file = resolve_input_manifest(paths_file, output, label="bins_file")

    if not os.path.isfile(paths_file):
            raise FileNotFoundError(f"Bin file not found: {paths_file}")

    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)
    errors = []

    with open(paths_file, "r") as f:
        for line in f:
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
                except DownloadError as exc:
                    errors.append(str(exc))
                    continue

            filename = fasta_re.sub("", os.path.basename(full_path))
            fasta_dict[filename] = full_path

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def path_bins_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"❌ Folder not found: {folder_path}")

    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)

    # Iterate over all files in the folder
    for file_name in os.listdir(folder_path):
        if fasta_re.search(file_name):
            full_path = os.path.join(folder_path, file_name)
            file_id = fasta_re.sub("", file_name)
            fasta_dict[file_id] = full_path

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

# Create dictionary of bin names and paths from the input file
def file_mags_to_json(paths_file=None, output=False):
    fasta_dict = {}
    fasta_re = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)

    paths_file = resolve_input_manifest(paths_file, output, label="mags_file")

    if not os.path.isfile(paths_file):
        raise FileNotFoundError(f"MAG file not found: {paths_file}")

    errors = []

    with open(paths_file, "r") as f:
        for line in f:
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
                except DownloadError as exc:
                    errors.append(str(exc))
                    continue

            filename = fasta_re.sub("", os.path.basename(full_path))
            fasta_dict[filename] = full_path

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/mags_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

#updated to account for compressed genomes
def path_mags_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"❌ Folder not found: {folder_path}")

    # Compile a regex that matches .fa/.fna/.fasta, optionally followed by .gz
    FASTA_RE = re.compile(r'\.(?:fa|fna|fasta)(?:\.gz)?$', re.IGNORECASE)

    # Iterate over all files in the folder
    fasta_dict = {}
    for fname in os.listdir(folder_path):
        if FASTA_RE.search(fname):
            full_path = os.path.join(folder_path, fname)
            sample_id = FASTA_RE.sub('', fname)
            fasta_dict[sample_id] = full_path

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
        except DownloadError as exc:
            errors.append(str(exc))

    if errors:
        print(f"ERROR: {len(errors)} file(s) could not be downloaded:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

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

    # Convert defaultdict to standard dict (optional)
    TRANSCRIPTOME_TO_READS1 = dict(TRANSCRIPTOME_TO_READS1)
    TRANSCRIPTOME_TO_READS2 = dict(TRANSCRIPTOME_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/transcriptome_to_reads1.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS1, f)

    with open(f"{output}/data/transcriptome_to_reads2.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS2, f)
