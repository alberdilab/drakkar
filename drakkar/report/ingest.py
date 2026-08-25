"""Loaders that project Drakkar output tables into the report database.

Every loader is independent and returns the number of rows it wrote, or None
when its source file is absent. Large tables are streamed in chunks so that a
catalogue-scale annotation table never has to fit in memory.
"""

import lzma
import math
from pathlib import Path

import pandas as pd
import yaml

from drakkar.report.schema import TAXONOMIC_RANKS, record_ingest

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


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def ingest_resources(connection, output_dir):
    """Load per-run metadata from the drakkar_<run_id>.yaml files."""
    output_path = Path(output_dir)
    metadata_files = sorted(output_path.glob("drakkar_*.yaml"))
    if not metadata_files:
        return None

    rows = []
    for path in metadata_files:
        try:
            with open(path, "r", encoding="utf-8") as handle:
                metadata = yaml.safe_load(handle) or {}
        except (OSError, yaml.YAMLError):
            continue
        if not isinstance(metadata, dict) or not metadata.get("run_id"):
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
    connection.commit()
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
        "singlem_fraction", "nonpareil_C", "nonpareil_LR",
        "nonpareil_modelR", "nonpareil_LRstar", "nonpareil_diversity",
    ]

    rows = []
    for row in frame.to_dict("records"):
        sample_id = _get(row, "sample")
        if sample_id is None:
            continue
        rows.append(
            (sample_id,)
            + tuple(_get(row, column, _int) for column in columns)
            + tuple(_get(row, column, _float) for column in float_columns)
        )

    placeholders = ", ".join(["?"] * (1 + len(columns) + len(float_columns)))
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
    connection.commit()
    return written


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
    connection.commit()
    return written


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
    connection.commit()
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
}
