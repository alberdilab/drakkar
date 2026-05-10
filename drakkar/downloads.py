import csv
import os
import shutil
import subprocess
import sys
import time
import re
from pathlib import Path
from urllib.parse import urlencode, urlparse
from urllib.request import urlopen

import pandas as pd

from drakkar.input_errors import (
    DownloadError,
    InputFileError,
    report_input_resolution_errors,
    require_non_empty_file,
)
from drakkar.output import print

REMOTE_URL_SCHEMES = {"http", "https", "ftp", "sftp"}

DEFAULT_DOWNLOAD_RETRIES = 5

READ1_BASENAME_PATTERN = re.compile(r"(?:^|[._-])(?:R?1)(?:[._-]|$)", re.IGNORECASE)

READ2_BASENAME_PATTERN = re.compile(r"(?:^|[._-])(?:R?2)(?:[._-]|$)", re.IGNORECASE)

def is_url(value):
    parsed = urlparse(str(value))
    return parsed.scheme in REMOTE_URL_SCHEMES

def _has_value(value):
    return not (pd.isna(value) or str(value).strip() == "")

def _normalized_value(value):
    if pd.isna(value):
        return ""
    return str(value).strip()

def _retry_delay(attempt):
    return 5 * (3 ** (attempt - 1))

def _download_urlopen(url, tmp_path):
    with urlopen(url) as response, open(tmp_path, "wb") as handle:
        shutil.copyfileobj(response, handle)

def _download_sftp(url, tmp_path):
    curl_path = shutil.which("curl")
    if not curl_path:
        raise DownloadError("sftp:// downloads require curl to be available in PATH.")
    result = subprocess.run(
        [curl_path, "--fail", "--location", "--silent", "--show-error", "--output", tmp_path, url],
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or f"curl exited with {result.returncode}"
        raise DownloadError(message)

def _download_once(url, tmp_path):
    parsed = urlparse(url)
    if parsed.scheme == "sftp":
        _download_sftp(url, tmp_path)
    else:
        _download_urlopen(url, tmp_path)


def _validate_paired_read_maps(reads1, reads2, input_dir, label):
    errors = []
    read1_samples = set(reads1)
    read2_samples = set(reads2)
    for sample in sorted(read2_samples - read1_samples):
        errors.append(f"{label} sample {sample} has reverse reads but no forward reads.")
    for sample in sorted(read1_samples - read2_samples):
        errors.append(f"{label} sample {sample} has forward reads but no reverse reads.")
    for sample in sorted(read1_samples & read2_samples):
        if not reads1[sample]:
            errors.append(f"{label} sample {sample} has no forward read files.")
        if not reads2[sample]:
            errors.append(f"{label} sample {sample} has no reverse read files.")
        for path in reads1[sample]:
            try:
                require_non_empty_file(path, f"{label} forward read for sample {sample}")
            except InputFileError as exc:
                errors.append(str(exc))
        for path in reads2[sample]:
            try:
                require_non_empty_file(path, f"{label} reverse read for sample {sample}")
            except InputFileError as exc:
                errors.append(str(exc))
    if not read1_samples and not read2_samples:
        errors.append(f"No paired read files matching *_1.fq.gz and *_2.fq.gz were found in {input_dir}.")
    if errors:
        report_input_resolution_errors(errors)


def download_to_cache(url, sample_name, column_name, output, cache_subdir="reads_cache", preserve_basename=False, max_retries=DEFAULT_DOWNLOAD_RETRIES):
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
        if os.path.getsize(dest_path) > 0:
            print(f"Using cached {column_name} for {sample_name}: {dest_path}", flush=True)
            return dest_path
        print(f"WARNING: Cached {column_name} for {sample_name} is empty and will be downloaded again: {dest_path}", flush=True)
        os.remove(dest_path)

    print(f"Downloading {column_name} for {sample_name} from {url}", flush=True)
    tmp_path = f"{dest_path}.tmp"
    for attempt in range(1, max_retries + 1):
        try:
            _download_once(url, tmp_path)
            if not os.path.isfile(tmp_path) or os.path.getsize(tmp_path) == 0:
                raise DownloadError("downloaded file is empty")
            os.replace(tmp_path, dest_path)
            print(f"Saved {column_name} for {sample_name} to {dest_path}", flush=True)
            return dest_path
        except Exception as exc:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            if attempt < max_retries:
                delay = _retry_delay(attempt)
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
    require_non_empty_file(resolved_path, f"{column_name} file on row {row_number}")
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
    for attempt in range(1, DEFAULT_DOWNLOAD_RETRIES + 1):
        try:
            with urlopen(query_url) as response:
                payload = response.read().decode("utf-8")
            if not payload.strip():
                raise DownloadError("ENA API returned an empty response")
            break
        except Exception as exc:
            if attempt < DEFAULT_DOWNLOAD_RETRIES:
                delay = _retry_delay(attempt)
                print(f"WARNING: ENA API query attempt {attempt}/{DEFAULT_DOWNLOAD_RETRIES} failed for accession {accession_value}: {exc}. Retrying in {delay}s...", flush=True)
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

def resolve_preprocessed_read_lists(row, row_number, output):
    """Like resolve_sample_read_lists but honours preprocessedreads1/2 columns first.

    Priority order:
      1. preprocessedreads1 / preprocessedreads2 columns (explicit preprocessed paths)
      2. preprocessing/final/<sample>_1.fq.gz (auto-detected from a prior preprocessing run)
      3. rawreads1 / rawreads2 or accession (raw inputs — same as resolve_sample_read_lists)
    """
    sample = row.get("sample")
    if not _has_value(sample):
        print(f"ERROR: Missing value in column 'sample' on row {row_number} of the sample info file.")
        sys.exit(1)

    sample_name = _normalized_value(sample)

    preprocessedreads1_value = _normalized_value(row.get("preprocessedreads1"))
    preprocessedreads2_value = _normalized_value(row.get("preprocessedreads2"))

    if preprocessedreads1_value or preprocessedreads2_value:
        for col, val in (("preprocessedreads1", preprocessedreads1_value), ("preprocessedreads2", preprocessedreads2_value)):
            if not val:
                print(f"ERROR: Column '{col}' is missing on row {row_number} while the paired column is present.")
                sys.exit(1)
        reads1 = [_resolve_reads_path(preprocessedreads1_value, sample_name, "preprocessedreads1", output, row_number)]
        reads2 = [_resolve_reads_path(preprocessedreads2_value, sample_name, "preprocessedreads2", output, row_number)]
        return sample_name, reads1, reads2

    final_r1 = os.path.join(output, "preprocessing", "final", f"{sample_name}_1.fq.gz")
    final_r2 = os.path.join(output, "preprocessing", "final", f"{sample_name}_2.fq.gz")
    if os.path.isfile(final_r1) and os.path.isfile(final_r2):
        require_non_empty_file(final_r1, f"auto-detected preprocessedreads1 for sample {sample_name}")
        require_non_empty_file(final_r2, f"auto-detected preprocessedreads2 for sample {sample_name}")
        return sample_name, [final_r1], [final_r2]

    return resolve_sample_read_lists(row, row_number, output)
