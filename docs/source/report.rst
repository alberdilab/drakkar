.. _report-database:

Report database reference
=========================

``drakkar report`` builds ``drakkar.db``, a SQLite projection of the tables
found in a Drakkar output directory. The database is the single source of truth
for the HTML report, so the ``.tsv`` and ``.tsv.xz`` outputs are read exactly
once, at ingest, and never again at render time.

The database is deliberately **not** a lossless copy. It is optimized for
querying and plotting; the source tables remain the archival record. In
particular, the ``details`` JSON of the annotation tables is not carried over.

Building the database
---------------------

.. code-block:: console

   $ drakkar report -o drakkar_output

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
- ``--primary-hits-only``: keep only rank-1 annotation hits. This roughly
  halves the largest table, at the cost of the secondary evidence that the 2.0
  annotation schema was designed to preserve. Leave it off if you intend to
  query ambiguous assignments.
- ``--list``: report which sections the output directory can support, without
  building anything.
- ``--force``: delete and rebuild an existing ``drakkar.db``.

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
