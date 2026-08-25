.. _report-database:

Report database and HTML report
===============================

``drakkar report`` does two things in one pass. It builds ``drakkar.db``, a
SQLite projection of the tables found in a Drakkar output directory, and it
renders ``drakkar_report.html`` from that database. The database is the single
source of truth for the report, so the ``.tsv`` and ``.tsv.xz`` outputs are
read exactly once, at ingest, and never again at render time.

The database is deliberately **not** a lossless copy. It is optimized for
querying and plotting; the source tables remain the archival record. In
particular, the ``details`` JSON of the annotation tables is not carried over.

Building the report
-------------------

.. code-block:: console

   $ drakkar report -o drakkar_output

This writes ``drakkar_output/drakkar.db`` and ``drakkar_output/drakkar_report.html``.

Because Drakkar is modular, a missing source file is the normal case rather
than an error. The command screens the output directory first and builds only
the sections whose inputs are present, naming the missing files for anything it
skips:

.. code-block:: console

   $ drakkar report -o drakkar_output --list

Options
-------

- ``-o/--output``: output directory to report on. Default is the current
  directory.
- ``--sections``: comma-separated list of sections to include. Valid values are
  ``preprocessing``, ``cataloging``, ``dereplication``, ``profiling``,
  ``taxonomy``, ``function``, ``expression``, ``resources``, and ``all``.
  Requesting a section whose inputs are absent warns and names the missing
  file rather than silently producing an empty section.
- ``--db-only``: build the database without rendering the HTML report.
- ``--html-only``: re-render ``drakkar_report.html`` from an existing
  ``drakkar.db`` without re-ingesting the source tables. Useful after
  upgrading Drakkar, or when the source directory is no longer at hand.
  Mutually exclusive with ``--db-only``.
- ``--primary-hits-only``: keep only rank-1 annotation hits. This roughly
  halves the largest table, at the cost of the secondary evidence that the 2.0
  annotation schema was designed to preserve. Leave it off if you intend to
  query ambiguous assignments.
- ``--list``: report which sections the output directory can support, without
  building anything.
- ``--force``: delete and rebuild an existing ``drakkar.db``. The HTML report
  is a derived artifact and is always overwritten, so it needs no flag.

The HTML report
---------------

``drakkar_report.html`` is a single self-contained file. The stylesheet and a
small navigation script are inlined and the Plotly bundle is embedded once, in
the first figure, so the report opens on a laptop with no network connection
and can be emailed or archived as-is. Nothing is fetched when the page loads.

The page is a left sidebar beside the sections it describes. The sidebar holds
the general information about the report — the Drakkar version that wrote the
database, the report schema version, the run identifiers found in the output
directory, the ingest timestamps, and which sections were rendered, which were
unavailable, and which were excluded by ``--sections`` — and, above it, the
table of contents. Sections that the database does not hold are named there
rather than silently dropped.

The table of contents is also the navigation: one section is shown at a time,
so the whole report is not laid out on a single scrolling page. Within a
section, tables longer than twenty rows are paged in the browser — every row is
still in the file, only twenty are on screen — and each table is preceded by
the averages of its numeric columns, shown as highlights. Printing, or opening
the page with scripting disabled, falls back to the flat page: every section
stacked and every row listed.

Each section is rendered from SQL aggregates only. ``gene_annotation``,
``cluster_annotation`` and ``gene_expression`` can hold tens of millions of
rows, so every figure and table on the page comes from a ``GROUP BY`` or a
``LIMIT``; no large table is ever read into memory whole. Long listings are
truncated with a note pointing back at the database, and the per-genome and
per-MAG heatmaps show the top entries by abundance or gene count.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Section
     - Contents
   * - Preprocessing
     - Per-sample read and base counts, metagenomic fraction, SingleM fraction
       and Nonpareil diversity, plus a stacked read-fate chart.
   * - Cataloging
     - Per-assembly QUAST and mapping statistics, bins contributed per binner,
       and per-sample mapping rates.
   * - Dereplication
     - Bin counts and mean quality before and after dereplication, and the ANI
       threshold used.
   * - Profiling
     - Genome quality summary, a completeness-against-contamination scatter,
       and a genome-by-sample relative abundance heatmap.
   * - Taxonomy
     - Distinct taxa and unclassified genomes at each GTDB rank, genomes per
       phylum, and phylum composition per sample when abundances are present.
   * - Functional annotation
     - Hits, annotated genes and distinct terms per source, a per-MAG
       annotation coverage heatmap, retained gene clusters and regions, and
       the annotation QC counts.
   * - Expression
     - Per-sample assigned counts and detected genes, and the length
       distribution of the quantified genes.
   * - Runs and resources
     - One row per recorded Drakkar run, with its command, modules, wall-clock
       duration and status.
   * - Provenance
     - The ingest log: every table, the file it came from, its row count and
       when it was ingested.

Rendering runs in-process in the plain Drakkar environment. It needs no
Snakemake and no Conda environment, so a report can be produced on a login node
or a laptop with only ``drakkar`` installed.

Sections and their inputs
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Section
     - Source files
   * - ``preprocessing``
     - ``preprocessing.tsv``
   * - ``cataloging``
     - ``cataloging.tsv``
   * - ``dereplication``
     - ``dereplicating.tsv``
   * - ``profiling``
     - ``profiling_genomes/final/{counts,bases,mags}.tsv``
   * - ``taxonomy``
     - ``annotating/genome_taxonomy.tsv``
   * - ``function``
     - ``annotating/gene_annotations.tsv.xz``,
       ``annotating/cluster_annotations.tsv.xz``,
       ``annotating/annotation_qc.tsv``
   * - ``expression``
     - ``expressing/gene_counts.tsv.xz``
   * - ``resources``
     - ``drakkar_<run_id>.yaml``

How the source tables are decluttered
-------------------------------------

Three transformations keep the database small enough to query comfortably at
catalogue scale.

**Repeated gene coordinates are lifted out.** ``gene_annotations.tsv.xz`` is a
long-form table, so ``mag``, ``gene``, ``contig``, ``start``, ``end`` and
``strand`` repeat on every hit row for the same gene. Those move into ``gene``,
which stores them once. ``source=prodigal`` rows are then omitted from
``gene_annotation`` entirely, because ``gene`` already records every predicted
gene.

**Repeated annotation labels are lifted out.** The human-readable ``annotation``
description repeats on every occurrence of an identifier — the same KO
description on every one of thousands of ``K00001`` rows. It moves into the
``annotation_term`` dimension table, keyed by ``(source, annotation_id)``.

**Packed fields are normalized.** The GTDB-Tk lineage string, the per-sample
mapping rates, the per-binner bin counts, and the wide genome-by-sample count
matrices each become their own table.

Taxonomy parsing
~~~~~~~~~~~~~~~~

``annotating/genome_taxonomy.tsv`` is a plain concatenation of the GTDB-Tk
bacterial and archaeal ``summary.tsv`` files, so it arrives with the lineage
packed into a single semicolon-joined ``classification`` string and several
free-text blob columns alongside it. Ingest splits the lineage into seven rank
columns, keeps the placement evidence, drops the blobs, and skips the archaeal
summary's repeated header row. GTDB-Tk renamed its closest-reference columns
between major versions, so both the ``fastani_*`` and ``closest_genome_*``
spellings are accepted.

Schema
------

Bookkeeping
~~~~~~~~~~~

``schema_version`` stamps the layout version and the Drakkar version that wrote
the database. ``ingest_log`` records, per table, the source file, its size and
modification time, the number of rows ingested, and when.

If the stored schema version does not match the running Drakkar build,
``drakkar report`` refuses to touch the database and asks for ``--force``.
Rebuilding from the source tables is cheap, so no migration path is provided.

Tables
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 26 48

   * - Table
     - Key
     - Contents
   * - ``run``
     - ``run_id``
     - One row per Drakkar invocation recorded in the output directory.
   * - ``sample``
     - ``sample_id``
     - Per-sample preprocessing statistics.
   * - ``assembly``
     - ``assembly_id``
     - Per-assembly QUAST, mapping and binning summary statistics.
   * - ``assembly_sample``
     - ``(assembly_id, sample_id)``
     - Which samples fed an assembly, which provided coverage, and the
       per-sample mapping rate.
   * - ``assembly_binner``
     - ``(assembly_id, binner)``
     - Bin count contributed by each binner.
   * - ``dereplication``
     - —
     - Single-row dereplication summary.
   * - ``genome``
     - ``genome_id``
     - Dereplicated genome metadata: quality, size, contiguity, cluster.
   * - ``genome_count``
     - ``(genome_id, sample_id)``
     - Long-form read counts and covered bases per genome and sample.
   * - ``genome_taxonomy``
     - ``genome_id``
     - Ranked GTDB lineage plus placement evidence.
   * - ``gene``
     - ``(mag, gene)``
     - One row per predicted gene, with its coordinates.
   * - ``annotation_term``
     - ``(source, annotation_id)``
     - One row per distinct annotation identifier and its label.
   * - ``gene_annotation``
     - ``(mag, gene, source, hit_rank)``
     - One row per accepted functional hit, excluding Prodigal rows.
   * - ``cluster_annotation``
     - ``(mag, source, cluster_id)``
     - One row per retained region or system.
   * - ``annotation_qc``
     - ``(mag, level, source)``
     - Per-MAG, per-source annotation QC counts.
   * - ``expressed_gene``
     - ``gene_id``
     - featureCounts gene identifiers and their coordinates.
   * - ``gene_expression``
     - ``(gene_id, sample_id)``
     - Long-form expression counts.

Expression identifiers come straight from featureCounts, which derives them
from the ``ID`` attribute of ``expressing/reference/metagenome.gff``. Drakkar
does not guarantee that they carry a MAG prefix, so the database does not
fabricate a join between ``gene_expression`` and ``gene``.

Querying the database
---------------------

The database is a plain SQLite file, so it can be queried from the shell, from
Python, or from R without Drakkar installed.

.. code-block:: console

   $ sqlite3 drakkar_output/drakkar.db '.tables'
   $ sqlite3 drakkar_output/drakkar.db 'SELECT phylum, COUNT(*) FROM genome_taxonomy GROUP BY phylum'

Genome abundance joined to taxonomy:

.. code-block:: sql

   SELECT t.phylum, c.sample_id, SUM(c.read_count) AS reads
   FROM genome_count AS c
   JOIN genome_taxonomy AS t ON t.genome_id = c.genome_id
   GROUP BY t.phylum, c.sample_id
   ORDER BY reads DESC;

The most frequent KEGG orthologs, with their descriptions:

.. code-block:: sql

   SELECT a.annotation_id, t.annotation, COUNT(*) AS genes
   FROM gene_annotation AS a
   JOIN annotation_term AS t
     ON t.source = a.source AND t.annotation_id = a.annotation_id
   WHERE a.source = 'kegg' AND a.is_primary = 1
   GROUP BY a.annotation_id, t.annotation
   ORDER BY genes DESC
   LIMIT 20;

From R:

.. code-block:: r

   library(RSQLite)
   con <- dbConnect(SQLite(), "drakkar_output/drakkar.db")
   counts <- dbGetQuery(con, "SELECT * FROM genome_count")
   taxonomy <- dbGetQuery(con, "SELECT * FROM genome_taxonomy")
   dbDisconnect(con)

From Python:

.. code-block:: python

   import sqlite3
   import pandas as pd

   with sqlite3.connect("drakkar_output/drakkar.db") as connection:
       counts = pd.read_sql("SELECT * FROM genome_count", connection)

For evidence that the database does not retain — the ``details`` JSON, native
alignment coordinates, source-specific scores — read the source tables
described in :doc:`annotation_tables` instead.
