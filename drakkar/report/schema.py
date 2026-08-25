"""SQLite schema for the Drakkar report database.

The database is a lean, queryable projection of the tabular outputs found in a
Drakkar output directory. It is deliberately not a lossless copy: verbose
provenance columns (notably the ``details`` JSON of the annotation tables) stay
in the source ``.tsv``/``.tsv.xz`` files, which remain the archival record.

Three decluttering rules shape the layout:

- Repeated gene coordinates are stored once in ``gene`` instead of once per
  annotation hit.
- Repeated annotation labels are stored once in ``annotation_term`` instead of
  once per hit.
- Packed multi-value fields (GTDB-Tk lineage strings, per-sample mapping rates,
  per-binner bin counts) are normalized into their own tables.
"""

import sqlite3
from datetime import datetime, timezone

# Bump when the layout below changes in a way that invalidates existing files.
SCHEMA_VERSION = 1

TAXONOMIC_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")

SCHEMA_STATEMENTS = [
    # -- Bookkeeping -------------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS schema_version (
        version INTEGER NOT NULL,
        drakkar_version TEXT,
        created_at TEXT NOT NULL
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS ingest_log (
        table_name TEXT PRIMARY KEY,
        section TEXT NOT NULL,
        source_file TEXT,
        source_mtime REAL,
        source_bytes INTEGER,
        rows_ingested INTEGER,
        ingested_at TEXT NOT NULL
    )
    """,
    # -- Run provenance ----------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS run (
        run_id TEXT PRIMARY KEY,
        drakkar_version TEXT,
        command TEXT,
        modules TEXT,
        started_at TEXT,
        finished_at TEXT,
        status TEXT,
        output_directory TEXT,
        argv TEXT
    )
    """,
    # -- Preprocessing -----------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS sample (
        sample_id TEXT PRIMARY KEY,
        seqkit_sana_reads INTEGER,
        reads_pre_fastp INTEGER,
        bases_pre_fastp INTEGER,
        adapter_trimmed_reads INTEGER,
        adapter_trimmed_bases INTEGER,
        reads_post_fastp INTEGER,
        bases_post_fastp INTEGER,
        host_reads INTEGER,
        host_bases INTEGER,
        metagenomic_reads INTEGER,
        metagenomic_bases INTEGER,
        singlem_fraction REAL,
        nonpareil_C REAL,
        nonpareil_LR REAL,
        nonpareil_modelR REAL,
        nonpareil_LRstar REAL,
        nonpareil_diversity REAL
    )
    """,
    # -- Cataloging --------------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS assembly (
        assembly_id TEXT PRIMARY KEY,
        assembly_contigs INTEGER,
        assembly_total_length INTEGER,
        assembly_largest_contig INTEGER,
        assembly_gc_percent REAL,
        assembly_N50 INTEGER,
        assembly_N75 INTEGER,
        assembly_L50 INTEGER,
        assembly_L75 INTEGER,
        mapped_reads INTEGER,
        total_reads INTEGER,
        mapping_rate_percent REAL,
        final_bins INTEGER,
        high_quality_bins INTEGER,
        medium_quality_bins INTEGER,
        low_quality_bins INTEGER,
        bin_total_size INTEGER,
        bin_mean_size REAL,
        bin_mean_N50 REAL,
        bin_total_contigs INTEGER,
        bin_mean_completeness REAL,
        bin_mean_contamination REAL,
        best_bin TEXT,
        best_bin_score REAL,
        best_bin_completeness REAL,
        best_bin_contamination REAL,
        best_bin_size INTEGER,
        best_bin_N50 INTEGER
    )
    """,
    # Normalized from the packed `samples` / `coverage_samples` columns.
    """
    CREATE TABLE IF NOT EXISTS assembly_sample (
        assembly_id TEXT NOT NULL,
        sample_id TEXT NOT NULL,
        is_assembly_input INTEGER NOT NULL DEFAULT 0,
        is_coverage_sample INTEGER NOT NULL DEFAULT 0,
        mapping_rate_percent REAL,
        PRIMARY KEY (assembly_id, sample_id)
    )
    """,
    # Normalized from the four `<binner>_bins` columns.
    """
    CREATE TABLE IF NOT EXISTS assembly_binner (
        assembly_id TEXT NOT NULL,
        binner TEXT NOT NULL,
        bin_count INTEGER,
        PRIMARY KEY (assembly_id, binner)
    )
    """,
    # -- Dereplication -----------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS dereplication (
        input_bin_number INTEGER,
        input_bin_completeness REAL,
        input_bin_contamination REAL,
        dereplication_ani REAL,
        output_bin_number INTEGER,
        output_bin_completeness REAL,
        output_bin_contamination REAL
    )
    """,
    # -- Profiling ---------------------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS genome (
        genome_id TEXT PRIMARY KEY,
        completeness REAL,
        contamination REAL,
        size INTEGER,
        contigs INTEGER,
        longest_contig INTEGER,
        n50 INTEGER,
        gc REAL,
        cluster TEXT,
        cluster_members INTEGER,
        score REAL
    )
    """,
    # Melted from the wide genome x sample counts.tsv / bases.tsv tables.
    """
    CREATE TABLE IF NOT EXISTS genome_count (
        genome_id TEXT NOT NULL,
        sample_id TEXT NOT NULL,
        read_count REAL,
        covered_bases REAL,
        PRIMARY KEY (genome_id, sample_id)
    )
    """,
    # -- Taxonomy ----------------------------------------------------------
    # Parsed from the raw GTDB-Tk summary concatenation: the semicolon-packed
    # `classification` string becomes seven rank columns, and the free-text
    # `other_related_references` / `note` blobs are dropped.
    """
    CREATE TABLE IF NOT EXISTS genome_taxonomy (
        genome_id TEXT PRIMARY KEY,
        domain TEXT,
        phylum TEXT,
        class TEXT,
        "order" TEXT,
        family TEXT,
        genus TEXT,
        species TEXT,
        classification_method TEXT,
        closest_reference TEXT,
        ani REAL,
        af REAL,
        red_value REAL,
        msa_percent REAL,
        warnings TEXT
    )
    """,
    # -- Functional annotation --------------------------------------------
    """
    CREATE TABLE IF NOT EXISTS gene (
        mag TEXT NOT NULL,
        gene TEXT NOT NULL,
        contig TEXT,
        start INTEGER,
        end INTEGER,
        strand TEXT,
        PRIMARY KEY (mag, gene)
    )
    """,
    # Dimension table: one row per distinct identifier instead of one per hit.
    """
    CREATE TABLE IF NOT EXISTS annotation_term (
        source TEXT NOT NULL,
        annotation_id TEXT NOT NULL,
        annotation TEXT,
        annotation_type TEXT,
        PRIMARY KEY (source, annotation_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS gene_annotation (
        mag TEXT NOT NULL,
        gene TEXT NOT NULL,
        source TEXT NOT NULL,
        hit_rank INTEGER NOT NULL,
        is_primary INTEGER,
        annotation_id TEXT,
        evalue REAL,
        bitscore REAL,
        rank_score REAL,
        PRIMARY KEY (mag, gene, source, hit_rank)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS cluster_annotation (
        mag TEXT NOT NULL,
        source TEXT NOT NULL,
        cluster_id TEXT NOT NULL,
        contig TEXT,
        start INTEGER,
        end INTEGER,
        type TEXT,
        gene_count INTEGER,
        substrate TEXT,
        pul_id TEXT,
        PRIMARY KEY (mag, source, cluster_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS annotation_qc (
        mag TEXT NOT NULL,
        level TEXT NOT NULL,
        source TEXT NOT NULL,
        reported_records INTEGER,
        retained_records INTEGER,
        rejected_records INTEGER,
        unmapped_records INTEGER,
        unique_entities INTEGER,
        filter_stage TEXT,
        PRIMARY KEY (mag, level, source)
    )
    """,
    # -- Expression --------------------------------------------------------
    # featureCounts identifiers are kept verbatim; Drakkar does not guarantee a
    # MAG prefix on them, so no MAG join is fabricated here.
    """
    CREATE TABLE IF NOT EXISTS expressed_gene (
        gene_id TEXT PRIMARY KEY,
        contig TEXT,
        start INTEGER,
        end INTEGER,
        strand TEXT,
        length INTEGER
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS gene_expression (
        gene_id TEXT NOT NULL,
        sample_id TEXT NOT NULL,
        count REAL,
        PRIMARY KEY (gene_id, sample_id)
    )
    """,
]

INDEX_STATEMENTS = [
    "CREATE INDEX IF NOT EXISTS idx_genome_count_sample ON genome_count (sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_genome_taxonomy_phylum ON genome_taxonomy (phylum)",
    "CREATE INDEX IF NOT EXISTS idx_genome_taxonomy_genus ON genome_taxonomy (genus)",
    "CREATE INDEX IF NOT EXISTS idx_gene_annotation_term ON gene_annotation (source, annotation_id)",
    "CREATE INDEX IF NOT EXISTS idx_gene_annotation_mag ON gene_annotation (mag)",
    "CREATE INDEX IF NOT EXISTS idx_cluster_annotation_source ON cluster_annotation (source)",
    "CREATE INDEX IF NOT EXISTS idx_gene_expression_sample ON gene_expression (sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_assembly_sample_sample ON assembly_sample (sample_id)",
]


def connect(db_path):
    """Open a SQLite connection tuned for bulk ingest."""
    connection = sqlite3.connect(str(db_path))
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode = WAL")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = MEMORY")
    connection.execute("PRAGMA cache_size = -262144")
    return connection


def create_schema(connection, drakkar_version=None):
    """Create every table and index, and stamp the schema version once."""
    for statement in SCHEMA_STATEMENTS:
        connection.execute(statement)
    for statement in INDEX_STATEMENTS:
        connection.execute(statement)
    existing = connection.execute("SELECT COUNT(*) FROM schema_version").fetchone()[0]
    if not existing:
        connection.execute(
            "INSERT INTO schema_version (version, drakkar_version, created_at) VALUES (?, ?, ?)",
            (SCHEMA_VERSION, drakkar_version, datetime.now(timezone.utc).isoformat()),
        )
    connection.commit()


def read_schema_version(connection):
    """Return the stored schema version, or None if the database predates it."""
    try:
        row = connection.execute(
            "SELECT version FROM schema_version ORDER BY rowid DESC LIMIT 1"
        ).fetchone()
    except sqlite3.DatabaseError:
        return None
    return row[0] if row else None


def record_ingest(connection, table_name, section, source_path, rows):
    """Record what was ingested into one table, for provenance and re-ingest."""
    mtime = None
    size = None
    if source_path is not None:
        try:
            stat = source_path.stat()
            mtime = stat.st_mtime
            size = stat.st_size
        except OSError:
            pass
    connection.execute(
        """
        INSERT INTO ingest_log
            (table_name, section, source_file, source_mtime, source_bytes, rows_ingested, ingested_at)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(table_name) DO UPDATE SET
            section = excluded.section,
            source_file = excluded.source_file,
            source_mtime = excluded.source_mtime,
            source_bytes = excluded.source_bytes,
            rows_ingested = excluded.rows_ingested,
            ingested_at = excluded.ingested_at
        """,
        (
            table_name,
            section,
            str(source_path) if source_path is not None else None,
            mtime,
            size,
            rows,
            datetime.now(timezone.utc).isoformat(),
        ),
    )


def finalize(connection):
    """Flush and compact the database after ingest."""
    connection.commit()
    connection.execute("PRAGMA optimize")
    connection.commit()
