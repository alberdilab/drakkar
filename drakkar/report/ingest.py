"""Loaders that project Drakkar output tables into the report database.

Every loader is independent and returns the number of rows it wrote, or None
when its source file is absent. Large tables are streamed in chunks so that a
catalogue-scale annotation table never has to fit in memory.
"""

import lzma
import math
import re
from pathlib import Path

import pandas as pd
import yaml

from drakkar.report.newick import count_tips
from drakkar.report.schema import TAXONOMIC_RANKS, record_ingest
from drakkar.report.sources import (
    benchmark_run_id,
    failure_run_id,
    find_benchmark_summaries,
    find_benchmark_tables,
    find_bin_reports,
    find_failure_reports,
    find_run_metadata,
)

# Rows per chunk when streaming the large annotation and expression tables.
CHUNK_SIZE = 200_000

# "-" is deliberately absent: it is a valid strand value, not a missing marker.
MISSING_TOKENS = {"", "NA", "N/A", "NAN", "NONE", "NULL"}

BINNER_COLUMNS = (
    ("metabat2_bins", "metabat2"),
    ("maxbin2_bins", "maxbin2"),
    ("semibin2_bins", "semibin2"),
    ("comebin_bins", "comebin"),
)

# Binette names each input bin set after the differing part of the paths it was
# given, which for Drakkar is the binner's rule directory: `metabat`, `maxbin`,
# `semibin`, `comebin`. Bins Binette assembled itself carry `binette` instead.
# The names are matched as substrings so that a single-binner run — where
# Binette falls back to the whole path — still resolves.
BINNER_ORIGIN_ALIASES = (
    ("metabat", "metabat2"),
    ("maxbin", "maxbin2"),
    ("semibin", "semibin2"),
    ("comebin", "comebin"),
    ("binette", "binette"),
)

# The name Binette gives bins it built itself, rather than took from a binner.
BINETTE_ORIGIN = "binette"

# GTDB-Tk renamed the closest-reference columns between major versions, so each
# logical field is looked up through a list of candidates.
GTDB_COLUMN_ALIASES = {
    "genome_id": ["user_genome", "genome", "magid"],
    "closest_reference": [
        "closest_genome_reference",
        "fastani_reference",
        "closest_placement_reference",
    ],
    "ani": ["closest_genome_ani", "fastani_ani", "closest_placement_ani"],
    "af": ["closest_genome_af", "fastani_af", "closest_placement_af"],
    "classification_method": ["classification_method"],
    "red_value": ["red_value"],
    "msa_percent": ["msa_percent"],
    "warnings": ["warnings"],
}

RANK_PREFIXES = {
    "d__": "domain",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}


# ---------------------------------------------------------------------------
# Value helpers
# ---------------------------------------------------------------------------

def _text(value):
    """Normalize a cell to a string, or None when it encodes a missing value."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    if text.upper() in MISSING_TOKENS:
        return None
    return text


def _float(value):
    text = _text(value)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _int(value):
    number = _float(value)
    if number is None or math.isnan(number) or math.isinf(number):
        return None
    return int(number)


def _fraction(value):
    """A share written either as ``0.85`` or as ``85%``, stored as a fraction.

    SingleM writes its read fraction with a percent sign in some versions and
    without one in others; the report reads the column as a fraction of one
    either way.
    """
    text = _text(value)
    if text is None:
        return None
    percent = text.endswith("%")
    number = _float(text[:-1] if percent else text)
    if number is None:
        return None
    return number / 100.0 if percent else number


def _bool(value):
    text = _text(value)
    if text is None:
        return None
    return 1 if text.lower() in {"true", "1", "yes", "t"} else 0


def _get(row, column, cast=_text):
    return cast(row[column]) if column in row else None


def _split_list(value):
    """Split a comma- or semicolon-packed field into its elements."""
    text = _text(value)
    if text is None:
        return []
    separator = ";" if ";" in text else ","
    return [item.strip() for item in text.split(separator) if item.strip()]


def _resolve(row, aliases):
    """Return the first present, non-missing value among candidate columns."""
    for column in aliases:
        if column in row:
            value = _text(row[column])
            if value is not None:
                return value
    return None


def _open_text(path):
    """Open a plain or xz-compressed text file."""
    if str(path).endswith(".xz"):
        return lzma.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _read_table(path, **kwargs):
    """Read a TSV, transparently handling xz compression."""
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("dtype", str)
    return pd.read_csv(path, **kwargs)


def _existing(output_dir, relative):
    path = Path(output_dir) / relative
    try:
        if path.is_file() and path.stat().st_size > 0:
            return path
    except OSError:
        pass
    return None


def _executemany(connection, statement, rows):
    if rows:
        connection.executemany(statement, rows)
    return len(rows)


def _read_csv(path, **kwargs):
    """Read one of dRep's comma-separated data tables."""
    kwargs.setdefault("sep", ",")
    return pd.read_csv(path, **kwargs)


def _genome_name(value):
    """Strip the directory and the FASTA extension off a dRep genome name.

    dRep writes bare file names in ``Cdb``/``Wdb`` but absolute paths in
    ``Ndb`` and ``Mdb``, so the two only join once both are reduced to the bin
    id Drakkar uses everywhere else.
    """
    text = _text(value)
    if text is None:
        return None
    name = Path(text).name
    for suffix in (".gz", ".fa", ".fna", ".fasta"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name or None


def _identity(value):
    """Read an ANI or a MASH similarity as a fraction of one.

    dRep writes these as fractions, but some versions of some comparison
    algorithms report percentages; anything above one is read as a percentage
    rather than stored as an impossible identity.
    """
    number = _float(value)
    if number is None:
        return None
    return number / 100.0 if number > 1.0 else number


def _canonical_binner(origin):
    """Map one Binette origin token onto the binner name the report uses."""
    text = _text(origin)
    if text is None:
        return None
    lowered = text.lower()
    for alias, canonical in BINNER_ORIGIN_ALIASES:
        if alias in lowered:
            return canonical
    return text


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

# Columns of the benchmark TSVs written by drakkar.benchmark, in file order,
# paired with the cast that turns a cell into a database value.
BENCHMARK_JOB_COLUMNS = (
    ("launch_index", _int),
    ("rule", _text),
    ("attempt", _int),
    ("logical_job_key", _text),
    ("internal_jobid", _text),
    ("external_jobid", _text),
    ("wildcards", _text),
    ("requested_cpus", _int),
    ("requested_mem_mb", _float),
    ("requested_runtime_min", _float),
    ("state", _text),
    ("exit_code", _text),
    ("alloc_cpus", _int),
    ("elapsed_sec", _float),
    ("cpu_time_sec", _float),
    ("max_rss_mb", _float),
    ("cpu_efficiency", _float),
    ("memory_efficiency", _float),
    ("runtime_efficiency", _float),
    ("oom", _bool),
    ("timeout", _bool),
)

BENCHMARK_RULE_COLUMNS = (
    ("rule", _text),
    ("launches", _int),
    ("logical_jobs", _int),
    ("retries", _int),
    ("failed_launches", _int),
    ("oom_launches", _int),
    ("timeout_launches", _int),
    ("median_requested_cpus", _float),
    ("median_alloc_cpus", _float),
    ("median_requested_mem_mb", _float),
    ("median_max_rss_mb", _float),
    ("median_memory_efficiency", _float),
    ("median_requested_runtime_min", _float),
    ("median_elapsed_sec", _float),
    ("median_runtime_efficiency", _float),
    ("allocated_cpu_sec", _float),
    ("used_cpu_sec", _float),
    ("weighted_cpu_efficiency", _float),
)

# Columns of the failure TSV drakkar.failures writes, in file order. The file
# also carries the run id in a column of its own, which is read separately so a
# table renamed away from its run still lands under the run that produced it.
FAILURE_COLUMNS = (
    ("rule", _text),
    ("target", _text),
    ("attempts", _int),
    ("status", _text),
    ("category", _text),
    ("reason", _text),
    ("slurm_state", _text),
    ("internal_jobid", _text),
    ("external_jobid", _text),
    ("detail", _text),
    ("action", _text),
    ("job_log", _text),
    ("output", _text),
    ("last_failure_at", _text),
)

# Run-level roll-up keys of drakkar_<run_id>.resources.yaml, after run_id.
BENCHMARK_SUMMARY_KEYS = (
    ("status", _text),
    ("profile", _text),
    ("message", _text),
    ("generated_at", _text),
    ("benchmarked_launches", _int),
    ("logical_jobs", _int),
    ("retries", _int),
    ("failed_launches", _int),
    ("oom_launches", _int),
    ("timeout_launches", _int),
    ("jobs_missing_accounting", _int),
    ("max_alloc_cpus", _int),
    ("peak_max_rss_mb", _float),
    ("total_elapsed_sec", _float),
    ("allocated_cpu_sec", _float),
    ("used_cpu_sec", _float),
    ("weighted_cpu_efficiency", _float),
)


def _load_yaml(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            document = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError):
        return None
    return document if isinstance(document, dict) else None


def ingest_resources(connection, output_dir):
    """Load per-run provenance and, when present, its resource benchmark.

    The provenance comes from the ``drakkar_<run_id>.yaml`` files. The usage
    figures come from the artefacts ``drakkar.benchmark`` writes for SLURM
    runs — a roll-up YAML per run and the per-launch / per-rule TSVs under
    ``benchmark/`` — and are simply absent for runs that were not benchmarked.
    The failure tables ``drakkar.failures`` writes are loaded alongside them,
    and are absent for a run in which nothing failed.
    """
    output_path = Path(output_dir)
    metadata_files = find_run_metadata(output_path)
    written = 0
    if metadata_files:
        written = _ingest_runs(connection, output_path, metadata_files)
    benchmark_rows = _ingest_benchmark(connection, output_path)
    failure_rows = _ingest_failures(connection, output_path)
    if not metadata_files and not benchmark_rows and not failure_rows:
        return None
    connection.commit()
    return written + benchmark_rows + failure_rows


def _ingest_runs(connection, output_path, metadata_files):
    """Write one `run` row per run metadata file."""
    rows = []
    for path in metadata_files:
        metadata = _load_yaml(path)
        if not metadata or not metadata.get("run_id"):
            continue
        modules = metadata.get("modules") or []
        argv = metadata.get("argv") or []
        rows.append(
            (
                _text(metadata.get("run_id")),
                _text(metadata.get("drakkar_version")),
                _text(metadata.get("command")),
                ",".join(str(item) for item in modules) if modules else None,
                _text(metadata.get("started_at")) or _text(metadata.get("timestamp")),
                _text(metadata.get("finished_at")),
                _text(metadata.get("status")),
                _text(metadata.get("output_directory")),
                " ".join(str(item) for item in argv) if argv else None,
            )
        )

    written = _executemany(
        connection,
        "INSERT OR REPLACE INTO run VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    record_ingest(connection, "run", "resources", output_path, written)
    return written


def _ingest_benchmark(connection, output_path):
    """Write the run-level, per-launch and per-rule benchmark rows.

    Returns the total number of rows written across the three tables, so that
    a directory holding benchmark files but no run metadata still counts as an
    ingested section.
    """
    written = 0
    written += _ingest_benchmark_summaries(connection, output_path)
    written += _ingest_benchmark_table(
        connection, output_path, "jobs", "benchmark_job", BENCHMARK_JOB_COLUMNS
    )
    written += _ingest_benchmark_table(
        connection, output_path, "rules", "benchmark_rule", BENCHMARK_RULE_COLUMNS
    )
    return written


def _ingest_benchmark_summaries(connection, output_path):
    summary_files = find_benchmark_summaries(output_path)
    if not summary_files:
        return 0

    rows = []
    for path in summary_files:
        summary = _load_yaml(path)
        if summary is None:
            continue
        run_id = _text(summary.get("run_id")) or benchmark_run_id(path)
        if run_id is None:
            continue
        rows.append(
            (run_id, *(cast(summary.get(key)) for key, cast in BENCHMARK_SUMMARY_KEYS))
        )

    placeholders = ", ".join(["?"] * (len(BENCHMARK_SUMMARY_KEYS) + 1))
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO run_benchmark VALUES ({placeholders})",
        rows,
    )
    record_ingest(connection, "run_benchmark", "resources", output_path, written)
    return written


def _ingest_benchmark_table(connection, output_path, kind, table, columns):
    """Load every ``logging/benchmark/drakkar_<run_id>.<kind>.tsv`` into one table."""
    sources = find_benchmark_tables(output_path, kind)
    if not sources:
        return 0

    rows = []
    for path in sources:
        run_id = benchmark_run_id(path, kind=kind)
        try:
            frame = _read_table(path)
        except (OSError, ValueError, pd.errors.ParserError):
            continue
        for record in frame.to_dict("records"):
            rows.append((run_id, *(_get(record, name, cast) for name, cast in columns)))

    placeholders = ", ".join(["?"] * (len(columns) + 1))
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO {table} VALUES ({placeholders})",
        rows,
    )
    record_ingest(connection, table, "resources", output_path / "benchmark", written)
    return written


def _ingest_failures(connection, output_path):
    """Load every ``logging/drakkar_<run_id>.failures.tsv`` into ``run_failure``.

    The rows are already one per rule and wildcard combination, carrying how
    many attempts that combination cost and whether a retry eventually
    recovered it, so they are stored as written rather than regrouped here.
    Their position in the file is kept as the key: two workflow-level errors of
    the same rule are distinct rows, and nothing else distinguishes them.
    """
    sources = find_failure_reports(output_path)
    if not sources:
        return 0

    rows = []
    for path in sources:
        run_id = failure_run_id(path)
        try:
            frame = _read_table(path)
        except (OSError, ValueError, pd.errors.ParserError):
            continue
        for index, record in enumerate(frame.to_dict("records"), start=1):
            rows.append(
                (
                    _get(record, "run_id") or run_id,
                    index,
                    *(_get(record, name, cast) for name, cast in FAILURE_COLUMNS),
                )
            )

    placeholders = ", ".join(["?"] * (len(FAILURE_COLUMNS) + 2))
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO run_failure VALUES ({placeholders})",
        rows,
    )
    record_ingest(connection, "run_failure", "resources", output_path, written)
    return written


def ingest_preprocessing(connection, output_dir):
    """Load preprocessing.tsv, one row per sample."""
    source = _existing(output_dir, "preprocessing.tsv")
    if source is None:
        return None

    frame = _read_table(source)
    columns = [
        "seqkit_sana_reads", "reads_pre_fastp", "bases_pre_fastp",
        "adapter_trimmed_reads", "adapter_trimmed_bases",
        "reads_post_fastp", "bases_post_fastp",
        "host_reads", "host_bases", "metagenomic_reads", "metagenomic_bases",
    ]
    float_columns = [
        "nonpareil_C", "nonpareil_LR", "nonpareil_modelR",
        "nonpareil_LRstar", "nonpareil_diversity",
    ]

    rows = []
    for row in frame.to_dict("records"):
        sample_id = _get(row, "sample")
        if sample_id is None:
            continue
        rows.append(
            (sample_id,)
            + tuple(_get(row, column, _int) for column in columns)
            + (_get(row, "singlem_fraction", _fraction),)
            + tuple(_get(row, column, _float) for column in float_columns)
        )

    placeholders = ", ".join(["?"] * (2 + len(columns) + len(float_columns)))
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO sample VALUES ({placeholders})",
        rows,
    )
    record_ingest(connection, "sample", "preprocessing", source, written)
    connection.commit()
    return written


def ingest_cataloging(connection, output_dir):
    """Load cataloging.tsv, normalizing its packed per-sample and binner fields."""
    source = _existing(output_dir, "cataloging.tsv")
    if source is None:
        return None

    frame = _read_table(source)
    int_columns = [
        "assembly_contigs", "assembly_total_length", "assembly_largest_contig",
    ]
    assembly_rows = []
    sample_rows = {}
    binner_rows = []

    for row in frame.to_dict("records"):
        assembly_id = _get(row, "assembly")
        if assembly_id is None:
            continue

        assembly_rows.append((
            assembly_id,
            _get(row, "assembly_contigs", _int),
            _get(row, "assembly_total_length", _int),
            _get(row, "assembly_largest_contig", _int),
            _get(row, "assembly_gc_percent", _float),
            _get(row, "assembly_N50", _int),
            _get(row, "assembly_N75", _int),
            _get(row, "assembly_L50", _int),
            _get(row, "assembly_L75", _int),
            _get(row, "mapped_reads", _int),
            _get(row, "total_reads", _int),
            _get(row, "mapping_rate_percent", _float),
            _get(row, "final_bins", _int),
            _get(row, "high_quality_bins", _int),
            _get(row, "medium_quality_bins", _int),
            _get(row, "low_quality_bins", _int),
            _get(row, "bin_total_size", _int),
            _get(row, "bin_mean_size", _float),
            _get(row, "bin_mean_N50", _float),
            _get(row, "bin_total_contigs", _int),
            _get(row, "bin_mean_completeness", _float),
            _get(row, "bin_mean_contamination", _float),
            _get(row, "best_bin"),
            _get(row, "best_bin_score", _float),
            _get(row, "best_bin_completeness", _float),
            _get(row, "best_bin_contamination", _float),
            _get(row, "best_bin_size", _int),
            _get(row, "best_bin_N50", _int),
        ))

        # `samples` and `coverage_samples` are packed lists; `sample_mapping_rates`
        # is packed as "sample:rate;sample:rate".
        for sample_id in _split_list(row.get("samples")):
            entry = sample_rows.setdefault((assembly_id, sample_id), [0, 0, None])
            entry[0] = 1
        for sample_id in _split_list(row.get("coverage_samples")):
            entry = sample_rows.setdefault((assembly_id, sample_id), [0, 0, None])
            entry[1] = 1
        for item in _split_list(row.get("sample_mapping_rates")):
            if ":" not in item:
                continue
            sample_id, _, rate = item.partition(":")
            sample_id = sample_id.strip()
            if not sample_id:
                continue
            entry = sample_rows.setdefault((assembly_id, sample_id), [0, 0, None])
            entry[2] = _float(rate)

        for column, binner in BINNER_COLUMNS:
            count = _get(row, column, _int)
            if count is not None:
                binner_rows.append((assembly_id, binner, count))

    placeholders = ", ".join(["?"] * 28)
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO assembly VALUES ({placeholders})",
        assembly_rows,
    )
    _executemany(
        connection,
        "INSERT OR REPLACE INTO assembly_sample VALUES (?, ?, ?, ?, ?)",
        [(a, s, v[0], v[1], v[2]) for (a, s), v in sample_rows.items()],
    )
    _executemany(
        connection,
        "INSERT OR REPLACE INTO assembly_binner VALUES (?, ?, ?)",
        binner_rows,
    )
    record_ingest(connection, "assembly", "cataloging", source, written)
    record_ingest(connection, "assembly_sample", "cataloging", source, len(sample_rows))
    record_ingest(connection, "assembly_binner", "cataloging", source, len(binner_rows))
    _ingest_bin_reports(connection, output_dir)
    connection.commit()
    return written


_BINETTE_BIN_NUMBER = re.compile(r"bin(\d+)$")


def _bin_name(assembly_id, reported):
    """Rebuild the Drakkar bin id from the name Binette gave the bin.

    Binette names its output `<prefix>_bin<n>`; Drakkar renames the FASTA to
    `<assembly>_bin_<n>.fa` and uses that id everywhere downstream, so the
    report stores the Drakkar form and keeps the two joinable.
    """
    text = _text(reported)
    if text is None:
        return None
    match = _BINETTE_BIN_NUMBER.search(text)
    return f"{assembly_id}_bin_{match.group(1)}" if match else text


def _ingest_bin_reports(connection, output_dir):
    """Load the per-assembly Binette reports behind ``cataloging/final``.

    They are what says where each final bin came from: `origin` names the
    binner that produced it, or `binette` when Binette assembled it out of
    contigs that no single binner had grouped that way. Without them the
    report can only say how many bins each binner made, not what became of
    them. Absent reports are the normal case for an output directory written
    before the reports were kept, and are not an error.
    """
    reports = find_bin_reports(output_dir)
    if not reports:
        return 0, 0

    bin_rows = []
    origin_rows = []
    for report in reports:
        assembly_id = report.stem
        try:
            frame = _read_table(report)
        except (OSError, ValueError, pd.errors.ParserError):
            continue
        for row in frame.to_dict("records"):
            bin_name = _bin_name(assembly_id, row.get("name") or row.get("bin_id"))
            if bin_name is None:
                continue
            origin = _get(row, "origin")
            is_original = _get(row, "is_original", _bool)
            bin_rows.append((
                assembly_id,
                bin_name,
                is_original,
                origin,
                _get(row, "original_name"),
                _get(row, "completeness", _float),
                _get(row, "contamination", _float),
                _get(row, "score", _float),
                _get(row, "size", _int),
                _get(row, "N50", _int),
                _get(row, "contig_count", _int),
            ))
            # Binette packs the origins of a bin several binners found
            # identically as "metabat;maxbin", so one bin can name two binners.
            for token in _split_list(origin):
                binner = _canonical_binner(token)
                if binner is not None:
                    origin_rows.append((assembly_id, bin_name, binner))

    written = _executemany(
        connection,
        "INSERT OR REPLACE INTO assembly_bin VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        bin_rows,
    )
    origins = _executemany(
        connection,
        "INSERT OR REPLACE INTO assembly_bin_origin VALUES (?, ?, ?)",
        sorted(set(origin_rows)),
    )
    if reports:
        record_ingest(connection, "assembly_bin", "cataloging", reports[0], written)
        record_ingest(connection, "assembly_bin_origin", "cataloging", reports[0], origins)
    return written, origins


def ingest_dereplication(connection, output_dir):
    """Load the single-row dereplicating.tsv summary."""
    source = _existing(output_dir, "dereplicating.tsv")
    if source is None:
        return None

    frame = _read_table(source)
    rows = [
        (
            _get(row, "input_bin_number", _int),
            _get(row, "input_bin_completeness", _float),
            _get(row, "input_bin_contamination", _float),
            _get(row, "dereplication_ani", _float),
            _get(row, "output_bin_number", _int),
            _get(row, "output_bin_completeness", _float),
            _get(row, "output_bin_contamination", _float),
        )
        for row in frame.to_dict("records")
    ]
    connection.execute("DELETE FROM dereplication")
    written = _executemany(
        connection,
        "INSERT INTO dereplication VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    record_ingest(connection, "dereplication", "dereplication", source, written)
    _ingest_drep_tables(connection, output_dir)
    connection.commit()
    return written


DREP_DATA_TABLES = "dereplicating/drep/data_tables"

# dRep writes bare file names in some tables and absolute paths in others, so
# every genome column is reduced to a bin id before anything is joined on it.
_FASTA_SUFFIX = re.compile(r"\.(fa|fna|fasta)$", re.IGNORECASE)
_GZIP_SUFFIX = re.compile(r"\.gz$", re.IGNORECASE)


def _genome_column(values):
    """Reduce a whole dRep genome column to bin ids at once.

    ``Mdb`` is quadratic in the number of bins, so this is done with vectorized
    string operations rather than a Python call per cell.
    """
    names = values.astype(str).str.rsplit("/", n=1).str[-1]
    names = names.str.replace(_GZIP_SUFFIX, "", regex=True)
    return names.str.replace(_FASTA_SUFFIX, "", regex=True)


def _nearest_mash_distances(path):
    """Return each bin's MASH distance to its closest other bin.

    Only the minimum is kept. The full matrix is quadratic and none of it is
    needed once the question is whether a bin had any near neighbour at all —
    which is what decides whether dereplication could have acted on it.
    """
    nearest = {}
    for chunk in _read_csv(path, chunksize=CHUNK_SIZE):
        columns = set(chunk.columns)
        if not {"genome1", "genome2"}.issubset(columns):
            return {}
        if "dist" in columns:
            distances = pd.to_numeric(chunk["dist"], errors="coerce")
        elif "similarity" in columns:
            distances = 1.0 - pd.to_numeric(chunk["similarity"], errors="coerce")
        else:
            return {}
        first = _genome_column(chunk["genome1"])
        second = _genome_column(chunk["genome2"])
        frame = pd.DataFrame({"genome": first, "other": second, "dist": distances})
        frame = frame[(frame["genome"] != frame["other"]) & frame["dist"].notna()]
        if frame.empty:
            continue
        for genome, distance in frame.groupby("genome")["dist"].min().items():
            value = float(distance)
            if genome not in nearest or value < nearest[genome]:
                nearest[genome] = value
    return nearest


def _pairwise_comparisons(path):
    """Return one averaged row per unordered pair of dRep's ANI comparisons.

    dRep compares each pair in both directions and clusters on the average of
    the two, so the report stores that average: it is the number the
    dereplication threshold was actually applied to, and storing both
    directions separately would double every bar of the histogram.
    """
    totals = {}
    for chunk in _read_csv(path, chunksize=CHUNK_SIZE):
        columns = set(chunk.columns)
        if not {"reference", "querry", "ani"}.issubset(columns):
            return []
        frame = pd.DataFrame({
            "a": _genome_column(chunk["reference"]),
            "b": _genome_column(chunk["querry"]),
            "ani": pd.to_numeric(chunk["ani"], errors="coerce"),
            "coverage": (
                pd.to_numeric(chunk["alignment_coverage"], errors="coerce")
                if "alignment_coverage" in columns
                else pd.Series([None] * len(chunk), index=chunk.index)
            ),
            "cluster": (
                chunk["primary_cluster"].astype(str)
                if "primary_cluster" in columns
                else pd.Series([None] * len(chunk), index=chunk.index)
            ),
        })
        frame = frame[(frame["a"] != frame["b"]) & frame["ani"].notna()]
        if frame.empty:
            continue
        # Order each pair so that A-B and B-A land on the same key, then fold
        # the chunk down in pandas: only the collapsed pairs cross into Python.
        frame["first"] = frame[["a", "b"]].min(axis=1)
        frame["second"] = frame[["a", "b"]].max(axis=1)
        grouped = frame.groupby(["first", "second"], sort=False).agg(
            ani_sum=("ani", "sum"),
            ani_count=("ani", "size"),
            coverage_sum=("coverage", "sum"),
            coverage_count=("coverage", "count"),
            cluster=("cluster", "first"),
        )
        for key, record in grouped.iterrows():
            entry = totals.setdefault(key, [0.0, 0, 0.0, 0, record["cluster"]])
            entry[0] += float(record["ani_sum"])
            entry[1] += int(record["ani_count"])
            entry[2] += float(record["coverage_sum"] or 0.0)
            entry[3] += int(record["coverage_count"])
    rows = []
    for (first, second), (ani_sum, ani_count, cov_sum, cov_count, cluster) in totals.items():
        ani = ani_sum / ani_count
        rows.append((
            first,
            second,
            _text(cluster),
            ani / 100.0 if ani > 1.0 else ani,
            (cov_sum / cov_count) if cov_count else None,
        ))
    return rows


def _ingest_drep_tables(connection, output_dir):
    """Load dRep's cluster assignments and pairwise comparisons.

    These are what turn "11 bins in, 6 out" into an account of why: which bins
    were near enough to anything else to be compared at all, how similar the
    compared pairs were, and where the ANI threshold fell among them. All four
    tables are optional — an output directory written by an older run, or one
    whose dRep working directory was cleaned up, still renders the summary.
    """
    cdb = _existing(output_dir, f"{DREP_DATA_TABLES}/Cdb.csv")
    if cdb is None:
        return 0, 0

    clusters = {}
    for row in _read_csv(cdb, dtype=str).to_dict("records"):
        genome = _genome_name(row.get("genome"))
        if genome is None:
            continue
        clusters[genome] = [
            _get(row, "primary_cluster"),
            _get(row, "secondary_cluster"),
            0,
            None,
        ]

    wdb = _existing(output_dir, f"{DREP_DATA_TABLES}/Wdb.csv")
    if wdb is not None:
        for row in _read_csv(wdb, dtype=str).to_dict("records"):
            genome = _genome_name(row.get("genome"))
            if genome in clusters:
                clusters[genome][2] = 1

    mdb = _existing(output_dir, f"{DREP_DATA_TABLES}/Mdb.csv")
    if mdb is not None:
        for genome, distance in _nearest_mash_distances(mdb).items():
            if genome in clusters:
                clusters[genome][3] = distance

    connection.execute("DELETE FROM genome_cluster")
    written = _executemany(
        connection,
        "INSERT INTO genome_cluster VALUES (?, ?, ?, ?, ?)",
        [(genome,) + tuple(values) for genome, values in sorted(clusters.items())],
    )
    record_ingest(connection, "genome_cluster", "dereplication", cdb, written)

    ndb = _existing(output_dir, f"{DREP_DATA_TABLES}/Ndb.csv")
    comparisons = 0
    if ndb is not None:
        connection.execute("DELETE FROM genome_comparison")
        comparisons = _executemany(
            connection,
            "INSERT INTO genome_comparison VALUES (?, ?, ?, ?, ?)",
            _pairwise_comparisons(ndb),
        )
        record_ingest(connection, "genome_comparison", "dereplication", ndb, comparisons)
    return written, comparisons


def ingest_profiling(connection, output_dir):
    """Load MAG metadata and melt the wide counts/bases tables into long form."""
    mags_source = _existing(output_dir, "profiling_genomes/final/mags.tsv")
    counts_source = _existing(output_dir, "profiling_genomes/final/counts.tsv")
    bases_source = _existing(output_dir, "profiling_genomes/final/bases.tsv")
    if mags_source is None and counts_source is None:
        return None

    written = 0
    if mags_source is not None:
        frame = _read_table(mags_source)
        rows = []
        for row in frame.to_dict("records"):
            genome_id = _get(row, "magid") or _get(row, "genome")
            if genome_id is None:
                continue
            rows.append((
                genome_id,
                _get(row, "completeness", _float),
                _get(row, "contamination", _float),
                _get(row, "size", _int),
                _get(row, "contigs", _int),
                _get(row, "longest_contig", _int),
                _get(row, "n50", _int),
                _get(row, "gc", _float),
                _get(row, "cluster"),
                _get(row, "cluster_members", _int),
                _get(row, "score", _float),
            ))
        written = _executemany(
            connection,
            "INSERT OR REPLACE INTO genome VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        record_ingest(connection, "genome", "profiling", mags_source, written)

    counts = _melt_matrix(counts_source)
    bases = _melt_matrix(bases_source)
    if counts or bases:
        keys = set(counts) | set(bases)
        rows = [
            (genome_id, sample_id, counts.get(key), bases.get(key))
            for key in sorted(keys)
            for genome_id, sample_id in (key,)
        ]
        _executemany(
            connection,
            "INSERT OR REPLACE INTO genome_count VALUES (?, ?, ?, ?)",
            rows,
        )
        record_ingest(
            connection, "genome_count", "profiling", counts_source or bases_source, len(rows)
        )

    connection.commit()
    return written


def _melt_matrix(source):
    """Melt a wide genome x sample matrix into a {(genome, sample): value} map."""
    if source is None:
        return {}
    frame = _read_table(source)
    if frame.empty or len(frame.columns) < 2:
        return {}
    genome_column = frame.columns[0]
    melted = {}
    for row in frame.to_dict("records"):
        genome_id = _text(row[genome_column])
        if genome_id is None:
            continue
        for column in frame.columns[1:]:
            value = _float(row[column])
            if value is not None:
                melted[(genome_id, str(column))] = value
    return melted


def ingest_taxonomy(connection, output_dir):
    """Parse the raw GTDB-Tk summary into ranked taxonomy columns.

    ``genome_taxonomy.tsv`` is a plain concatenation of the GTDB-Tk bacterial
    and archaeal summaries, so the lineage arrives as one semicolon-packed
    ``classification`` string alongside several free-text blobs. Only the
    ranked lineage and the placement evidence are kept.
    """
    source = _existing(output_dir, "annotating/genome_taxonomy.tsv")
    if source is None:
        return None

    frame = _read_table(source)
    rows = []
    for row in frame.to_dict("records"):
        genome_id = _resolve(row, GTDB_COLUMN_ALIASES["genome_id"])
        if genome_id is None:
            continue
        # Concatenating the two GTDB-Tk summaries repeats the header row.
        if genome_id.lower() == "user_genome":
            continue
        ranks = _parse_lineage(row.get("classification"))
        rows.append(
            (genome_id,)
            + tuple(ranks.get(rank) for rank in TAXONOMIC_RANKS)
            + (
                _resolve(row, GTDB_COLUMN_ALIASES["classification_method"]),
                _resolve(row, GTDB_COLUMN_ALIASES["closest_reference"]),
                _float(_resolve(row, GTDB_COLUMN_ALIASES["ani"])),
                _float(_resolve(row, GTDB_COLUMN_ALIASES["af"])),
                _float(_resolve(row, GTDB_COLUMN_ALIASES["red_value"])),
                _float(_resolve(row, GTDB_COLUMN_ALIASES["msa_percent"])),
                _resolve(row, GTDB_COLUMN_ALIASES["warnings"]),
            )
        )

    placeholders = ", ".join(["?"] * (1 + len(TAXONOMIC_RANKS) + 7))
    written = _executemany(
        connection,
        f"INSERT OR REPLACE INTO genome_taxonomy VALUES ({placeholders})",
        rows,
    )
    record_ingest(connection, "genome_taxonomy", "taxonomy", source, written)
    _ingest_taxonomy_trees(connection, output_dir)
    connection.commit()
    return written


# The pruned GTDB-Tk placement trees, and the domain each one speaks for. Both
# are optional: a catalogue with no archaea has no archaeal tree, and a run
# whose GTDB-Tk output predates the pruning rule has neither.
TAXONOMY_TREES = (
    ("bacteria", "annotating/bacteria.tree"),
    ("archaea", "annotating/archaea.tree"),
)

# A tree of one tip has no topology to draw and is stored as if it were absent.
MIN_TREE_TIPS = 2


def _ingest_taxonomy_trees(connection, output_dir):
    """Store the pruned placement trees as the Newick text they arrive as.

    The pruning rule writes an empty ``bacteria.tree`` when GTDB-Tk produced no
    classify tree at all, so a file is kept only once it parses into a tree
    with tips to lay out.
    """
    written = 0
    for domain, relative in TAXONOMY_TREES:
        source = _existing(output_dir, relative)
        if source is None:
            continue
        try:
            newick = source.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        tips = count_tips(newick)
        if tips < MIN_TREE_TIPS:
            continue
        connection.execute(
            "INSERT OR REPLACE INTO genome_tree VALUES (?, ?, ?)",
            (domain, newick, tips),
        )
        record_ingest(connection, f"genome_tree:{domain}", "taxonomy", source, 1)
        written += 1
    return written


def _parse_lineage(value):
    """Split a `d__X;p__Y;...` GTDB lineage into a rank -> name mapping."""
    ranks = {}
    text = _text(value)
    if text is None:
        return ranks
    for element in text.split(";"):
        element = element.strip()
        prefix = element[:3]
        rank = RANK_PREFIXES.get(prefix)
        if rank is None:
            continue
        name = element[3:].strip()
        ranks[rank] = name or None
    return ranks


def ingest_function(connection, output_dir, primary_hits_only=False):
    """Stream the annotation tables into their decluttered relational form.

    The long-form gene table repeats gene coordinates and annotation labels on
    every hit row. Those are lifted into ``gene`` and ``annotation_term``, the
    ``details`` JSON is dropped, and ``source=prodigal`` rows are omitted from
    ``gene_annotation`` because the ``gene`` table already records every
    predicted gene.
    """
    gene_source = _existing(output_dir, "annotating/gene_annotations.tsv.xz")
    cluster_source = _existing(output_dir, "annotating/cluster_annotations.tsv.xz")
    qc_source = _existing(output_dir, "annotating/annotation_qc.tsv")
    if gene_source is None and cluster_source is None:
        return None

    gene_rows_written = 0
    annotation_rows_written = 0
    term_rows_written = 0

    if gene_source is not None:
        seen_terms = set()
        with _open_text(gene_source) as handle:
            reader = pd.read_csv(handle, sep="\t", dtype=str, chunksize=CHUNK_SIZE)
            for chunk in reader:
                genes = {}
                annotations = []
                terms = {}
                for row in chunk.to_dict("records"):
                    mag = _get(row, "mag")
                    gene = _get(row, "gene")
                    if mag is None or gene is None:
                        continue
                    genes.setdefault((mag, gene), (
                        mag, gene,
                        _get(row, "contig"),
                        _get(row, "start", _int),
                        _get(row, "end", _int),
                        _get(row, "strand"),
                    ))

                    source_name = _get(row, "source")
                    # Prodigal rows carry no functional evidence; the gene table
                    # already records that the gene exists.
                    if source_name is None or source_name == "prodigal":
                        continue

                    is_primary = _get(row, "is_primary", _bool)
                    if primary_hits_only and is_primary != 1:
                        continue

                    annotation_id = _get(row, "annotation_id")
                    if annotation_id is not None:
                        key = (source_name, annotation_id)
                        if key not in seen_terms and key not in terms:
                            terms[key] = (
                                source_name,
                                annotation_id,
                                _get(row, "annotation"),
                                _get(row, "annotation_type"),
                            )

                    annotations.append((
                        mag, gene, source_name,
                        _get(row, "hit_rank", _int),
                        is_primary,
                        annotation_id,
                        _get(row, "evalue", _float),
                        _get(row, "bitscore", _float),
                        _get(row, "rank_score", _float),
                    ))

                gene_rows_written += _executemany(
                    connection,
                    "INSERT OR IGNORE INTO gene VALUES (?, ?, ?, ?, ?, ?)",
                    list(genes.values()),
                )
                term_rows_written += _executemany(
                    connection,
                    "INSERT OR IGNORE INTO annotation_term VALUES (?, ?, ?, ?)",
                    list(terms.values()),
                )
                seen_terms.update(terms)
                annotation_rows_written += _executemany(
                    connection,
                    "INSERT OR REPLACE INTO gene_annotation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    annotations,
                )
                connection.commit()

        record_ingest(connection, "gene", "function", gene_source, gene_rows_written)
        record_ingest(connection, "annotation_term", "function", gene_source, term_rows_written)
        record_ingest(connection, "gene_annotation", "function", gene_source, annotation_rows_written)

    if cluster_source is not None:
        with _open_text(cluster_source) as handle:
            frame = pd.read_csv(handle, sep="\t", dtype=str)
        rows = []
        for row in frame.to_dict("records"):
            mag = _get(row, "mag")
            source_name = _get(row, "source")
            cluster_id = _get(row, "cluster_id")
            if mag is None or source_name is None or cluster_id is None:
                continue
            rows.append((
                mag, source_name, cluster_id,
                _get(row, "contig"),
                _get(row, "start", _int),
                _get(row, "end", _int),
                _get(row, "type"),
                _get(row, "gene_count", _int),
                _get(row, "substrate"),
                _get(row, "pul_id"),
            ))
        written = _executemany(
            connection,
            "INSERT OR REPLACE INTO cluster_annotation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        record_ingest(connection, "cluster_annotation", "function", cluster_source, written)

    if qc_source is not None:
        frame = _read_table(qc_source)
        rows = []
        for row in frame.to_dict("records"):
            mag = _get(row, "mag")
            level = _get(row, "level")
            source_name = _get(row, "source")
            if mag is None or level is None or source_name is None:
                continue
            rows.append((
                mag, level, source_name,
                _get(row, "reported_records", _int),
                _get(row, "retained_records", _int),
                _get(row, "rejected_records", _int),
                _get(row, "unmapped_records", _int),
                _get(row, "unique_entities", _int),
                _get(row, "filter_stage"),
            ))
        written = _executemany(
            connection,
            "INSERT OR REPLACE INTO annotation_qc VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        record_ingest(connection, "annotation_qc", "function", qc_source, written)

    connection.commit()
    return annotation_rows_written


def ingest_expression(connection, output_dir):
    """Melt the featureCounts gene x sample matrix into long form."""
    source = _existing(output_dir, "expressing/gene_counts.tsv.xz")
    if source is None:
        return None

    skiprows = _count_leading_comments(source)
    fixed = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    gene_rows_written = 0
    count_rows_written = 0

    with _open_text(source) as handle:
        reader = pd.read_csv(
            handle, sep="\t", dtype=str, skiprows=skiprows, chunksize=CHUNK_SIZE
        )
        for chunk in reader:
            sample_columns = [
                column for column in chunk.columns if column not in fixed
            ]
            genes = []
            counts = []
            for row in chunk.to_dict("records"):
                gene_id = _get(row, "Geneid")
                if gene_id is None:
                    continue
                genes.append((
                    gene_id,
                    _first_element(row.get("Chr")),
                    _bound(row.get("Start"), min),
                    _bound(row.get("End"), max),
                    _first_element(row.get("Strand")),
                    _get(row, "Length", _int),
                ))
                for sample_id in sample_columns:
                    value = _float(row[sample_id])
                    if value is not None:
                        counts.append((gene_id, str(sample_id), value))

            gene_rows_written += _executemany(
                connection,
                "INSERT OR IGNORE INTO expressed_gene VALUES (?, ?, ?, ?, ?, ?)",
                genes,
            )
            count_rows_written += _executemany(
                connection,
                "INSERT OR REPLACE INTO gene_expression VALUES (?, ?, ?)",
                counts,
            )
            connection.commit()

    record_ingest(connection, "expressed_gene", "expression", source, gene_rows_written)
    record_ingest(connection, "gene_expression", "expression", source, count_rows_written)
    connection.commit()
    return count_rows_written


def _count_leading_comments(path):
    """Count the leading '#' lines that featureCounts writes before the header."""
    count = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                count += 1
            else:
                break
    return count


def _first_element(value):
    """featureCounts joins multi-exon attributes with ';'; take the first."""
    elements = _split_list(value)
    return elements[0] if elements else None


def _bound(value, selector):
    """Reduce a ';'-joined coordinate list to its min or max."""
    numbers = [_int(item) for item in _split_list(value)]
    numbers = [number for number in numbers if number is not None]
    return selector(numbers) if numbers else None


# ---------------------------------------------------------------------------
# Antimicrobial resistance
# ---------------------------------------------------------------------------

def ingest_amr(connection, output_dir):
    """Load the assembly-level AMR tables written by `aggregate_amr`.

    The per-assembly counts arrive as two files keyed on `assembly_id` — the
    summary and the QC roll-up — and are folded into one row, because every QC
    figure is the denominator of a summary figure. The larger per-hit and
    per-locus tables are streamed, and the verbose `raw_details` / `details`
    JSON columns are left in the source files.
    """
    summary_source = _existing(output_dir, "amr/assembly_summary.tsv")
    qc_source = _existing(output_dir, "amr/amr_qc.tsv")
    if summary_source is None and qc_source is None:
        return None

    written = _ingest_amr_assemblies(connection, summary_source, qc_source)
    written += _ingest_amr_hits(connection, output_dir)
    written += _ingest_amr_loci(connection, output_dir)
    written += _ingest_amr_drug_classes(connection, output_dir)
    written += _ingest_amr_mobility(connection, output_dir)
    connection.commit()
    return written


def _ingest_amr_assemblies(connection, summary_source, qc_source):
    """Merge assembly_summary.tsv and amr_qc.tsv into one row per assembly."""
    merged = {}
    for source, columns in (
        (summary_source, ("assembly_type", "organism", "contig_count",
                          "total_length", "amrfinder_hits", "rgi_hits",
                          "mobility_regions", "amr_loci", "multi_tool_loci",
                          "mobile_loci")),
        (qc_source, ("amrfinder_hits", "rgi_hits", "mobility_regions",
                     "amr_loci", "multi_tool_loci", "mobile_loci",
                     "amrfinder_hits_without_coordinates",
                     "rgi_hits_without_coordinates", "mobility_links")),
    ):
        if source is None:
            continue
        for row in _read_table(source).to_dict("records"):
            assembly_id = _get(row, "assembly_id")
            if assembly_id is None:
                continue
            record = merged.setdefault(assembly_id, {})
            for column in columns:
                value = _get(row, column, _text if column in {"assembly_type", "organism"} else _int)
                if value is not None:
                    record[column] = value

    rows = [
        (
            assembly_id,
            record.get("assembly_type"),
            record.get("organism"),
            record.get("contig_count"),
            record.get("total_length"),
            record.get("amrfinder_hits"),
            record.get("rgi_hits"),
            record.get("mobility_regions"),
            record.get("amr_loci"),
            record.get("multi_tool_loci"),
            record.get("mobile_loci"),
            record.get("amrfinder_hits_without_coordinates"),
            record.get("rgi_hits_without_coordinates"),
            record.get("mobility_links"),
        )
        for assembly_id, record in sorted(merged.items())
    ]
    written = _executemany(
        connection,
        "INSERT OR REPLACE INTO amr_assembly VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    record_ingest(
        connection, "amr_assembly", "amr", summary_source or qc_source, written
    )
    return written


def _stream_amr_table(connection, output_dir, relative, table, statement, build):
    """Stream one AMR table into the database, `build` mapping row -> tuple."""
    source = _existing(output_dir, relative)
    if source is None:
        return 0
    written = 0
    with _open_text(source) as handle:
        for chunk in pd.read_csv(handle, sep="\t", dtype=str, chunksize=CHUNK_SIZE):
            rows = [
                tuple_ for tuple_ in (build(row) for row in chunk.to_dict("records"))
                if tuple_ is not None
            ]
            written += _executemany(connection, statement, rows)
            connection.commit()
    record_ingest(connection, table, "amr", source, written)
    return written


def _ingest_amr_hits(connection, output_dir):
    def build(row):
        assembly_id = _get(row, "assembly_id")
        hit_id = _get(row, "hit_id")
        if assembly_id is None or hit_id is None:
            return None
        return (
            assembly_id,
            hit_id,
            _get(row, "locus_id"),
            _get(row, "source"),
            _get(row, "contig"),
            _get(row, "start", _int),
            _get(row, "end", _int),
            _get(row, "strand"),
            _get(row, "gene_symbol"),
            _get(row, "gene_name"),
            _get(row, "ontology_id"),
            _get(row, "drug_class"),
            _get(row, "resistance_mechanism"),
            _get(row, "gene_family"),
            _get(row, "detection_grade"),
            _get(row, "method"),
            _get(row, "identity", _float),
            _get(row, "reference_coverage", _float),
            _get(row, "bitscore", _float),
            _get(row, "is_partial", _bool),
        )

    return _stream_amr_table(
        connection, output_dir, "amr/amr_hits.tsv.xz", "amr_hit",
        "INSERT OR REPLACE INTO amr_hit VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        build,
    )


def _ingest_amr_loci(connection, output_dir):
    def build(row):
        assembly_id = _get(row, "assembly_id")
        locus_id = _get(row, "locus_id")
        if assembly_id is None or locus_id is None:
            return None
        return (
            assembly_id,
            locus_id,
            _get(row, "contig"),
            _get(row, "start", _int),
            _get(row, "end", _int),
            _get(row, "strand"),
            _get(row, "primary_gene"),
            _get(row, "gene_symbols"),
            _get(row, "gene_families"),
            _get(row, "ontology_ids"),
            _get(row, "drug_classes"),
            _get(row, "drug_subclasses"),
            _get(row, "resistance_mechanisms"),
            _get(row, "sources"),
            _get(row, "source_count", _int),
            _get(row, "hit_count", _int),
            _get(row, "support_status"),
            _get(row, "concordance"),
        )

    return _stream_amr_table(
        connection, output_dir, "amr/amr_loci.tsv.xz", "amr_locus",
        "INSERT OR REPLACE INTO amr_locus VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        build,
    )


def _ingest_amr_drug_classes(connection, output_dir):
    def build(row):
        assembly_id = _get(row, "assembly_id")
        hit_id = _get(row, "hit_id")
        # A hit whose caller named no drug class still gets a row in the source
        # table, with an empty class. It says nothing this table is for, and
        # its mechanism and gene family are already on `amr_hit`, so it is left
        # out rather than stored under a null class.
        drug_class = _get(row, "drug_class")
        if assembly_id is None or hit_id is None or drug_class is None:
            return None
        return (
            assembly_id,
            hit_id,
            drug_class,
            _get(row, "locus_id"),
            _get(row, "source"),
            _get(row, "drug_subclass"),
            _get(row, "resistance_mechanism"),
            _get(row, "gene_family"),
            _get(row, "ontology_id"),
        )

    return _stream_amr_table(
        connection, output_dir, "amr/amr_drug_classes.tsv.xz", "amr_drug_class",
        "INSERT OR REPLACE INTO amr_drug_class VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        build,
    )


def _ingest_amr_mobility(connection, output_dir):
    def build_region(row):
        assembly_id = _get(row, "assembly_id")
        region_id = _get(row, "region_id")
        if assembly_id is None or region_id is None:
            return None
        return (
            assembly_id,
            region_id,
            _get(row, "contig"),
            _get(row, "start", _int),
            _get(row, "end", _int),
            _get(row, "context_type"),
            _get(row, "seq_name"),
            _get(row, "length", _int),
            _get(row, "topology"),
            _get(row, "score", _float),
            _get(row, "fdr", _float),
            _get(row, "marker_enrichment", _float),
            _get(row, "hallmark_count", _int),
            _get(row, "gene_count", _int),
            _get(row, "conjugation_genes"),
            _get(row, "taxonomy"),
        )

    def build_link(row):
        assembly_id = _get(row, "assembly_id")
        locus_id = _get(row, "locus_id")
        region_id = _get(row, "region_id")
        if assembly_id is None or locus_id is None or region_id is None:
            return None
        return (
            assembly_id,
            locus_id,
            region_id,
            _get(row, "context_type"),
            _get(row, "contig"),
            _get(row, "overlap_bp", _int),
            _get(row, "locus_overlap_fraction", _float),
            _get(row, "region_score", _float),
        )

    written = _stream_amr_table(
        connection, output_dir, "amr/mobility_regions.tsv.xz", "amr_mobility_region",
        "INSERT OR REPLACE INTO amr_mobility_region VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        build_region,
    )
    written += _stream_amr_table(
        connection, output_dir, "amr/amr_mobility.tsv.xz", "amr_mobility",
        "INSERT OR REPLACE INTO amr_mobility VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        build_link,
    )
    return written


# Section name -> loader. `resources` is handled first so that run provenance
# is present even when every data section is missing.
SECTION_LOADERS = {
    "resources": ingest_resources,
    "preprocessing": ingest_preprocessing,
    "cataloging": ingest_cataloging,
    "dereplication": ingest_dereplication,
    "profiling": ingest_profiling,
    "taxonomy": ingest_taxonomy,
    "function": ingest_function,
    "expression": ingest_expression,
    "amr": ingest_amr,
}
