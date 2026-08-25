Workflow guide
==============

This page covers the analysis workflows in DRAKKAR: the complete pipeline and
the module-level commands used to run specific stages independently.

Workflow overview
-----------------

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Command
     - Purpose
     - Typical outputs
   * - ``drakkar complete``
     - Run the main workflow end-to-end from reads to downstream products.
     - Full output tree across preprocessing, cataloging, profiling,
       annotation, and optional expression.
   * - ``drakkar preprocessing``
     - Clean reads and optionally remove host DNA.
     - Cleaned read files, preprocessing summaries, microbial fraction, and
       optional Nonpareil outputs.
   * - ``drakkar cataloging``
     - Assemble reads, bin contigs, and summarize the MAG catalog.
     - Assemblies, bins, bin metadata, and ``cataloging.tsv``.
   * - ``drakkar profiling``
     - Dereplicate MAGs and quantify genomes or pangenomes across samples.
     - Dereplicated genomes, abundance tables, and profiling outputs.
   * - ``drakkar annotating``
     - Annotate MAGs taxonomically and functionally.
     - Taxonomy tables plus gene- and cluster-level annotation tables.
   * - ``drakkar expressing``
     - Map metatranscriptomes to annotated genes.
     - Gene expression tables under ``expressing/``.
   * - ``drakkar dereplicating``
     - Run dereplication only, without read mapping.
     - Dereplicated genomes in ``dereplicating/final``.
   * - ``drakkar inspecting``
     - Run microdiversity and mapping inspection steps.
     - Inspection outputs derived from MAGs, coverage tables, and BAM files.

Complete workflow
-----------------

Run the full pipeline in sequence:

.. code-block:: console

   $ drakkar complete -f input_info.tsv -o drakkar_output -m individual -t genomes

Options:

- ``-i/--input``: input directory for reads.
- ``-f/--file``: sample info table (TSV, CSV, or semicolon-separated), with
  read pairs provided either as ``rawreads1``/``rawreads2`` or as an ENA/SRA
  ``accession``.
- ``-o/--output``: output directory.
- ``-r/--reference``: local path, URL, or NCBI genome assembly accession such
  as ``GCF_000001405.40`` for a host reference genome used in preprocessing.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``-m/--mode``: assembly modes such as ``individual`` and ``all``.
- ``-b/--binners``: comma-separated binners for cataloging
  (``metabat``, ``maxbin``, ``semibin``, ``comebin``; default: all).
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``--annotation-type``: comma-separated annotation targets. See *Annotating*
  below for the full set.
- ``--annotation-evalue``: maximum fallback e-value for applicable merged gene
  annotation hits (default: ``1e-10``); sources with native cutoffs retain
  those policies.
- ``--annotation-identity``: minimum percent identity for merged gene
  annotation hits with identity values, currently VFDB/MMseqs hits
  (default: ``50``).
- ``--annotation-query-coverage``: minimum VFDB/MMseqs query coverage as a
  fraction from 0 to 1 (default: ``0.5``).
- ``--annotation-target-coverage``: minimum VFDB/MMseqs target coverage as a
  fraction from 0 to 1 (default: ``0.5``).
- ``--min-completeness`` / ``--max-contamination`` / ``--min-bin-length`` /
  ``--max-bin-length``: bin filtering thresholds applied by Binette during
  cataloging. See *Cataloging* below.
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``--fraction``: compute microbial fraction with SingleM.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil.
- ``-a/--ani``: dRep ANI threshold (default: ``0.98``).
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile (default: ``slurm``).
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark``: skip SLURM resource benchmark collection after the run.
- ``--memory-multiplier N`` / ``--time-multiplier N``: scale per-rule resource
  requests before the configured caps are applied.
- ``--snakemake-*`` / ``--slurm-*``: Snakemake and SLURM override flags. See
  :ref:`snakemake-slurm-management` in :doc:`operations`.

Module reference
----------------

Preprocessing
^^^^^^^^^^^^^

Quality filters reads, optionally removes host DNA, and writes cleaned reads
and preprocessing summaries.

.. code-block:: console

   $ drakkar preprocessing -i /path/to/reads -o drakkar_output -r host.fna

Options:

- ``-i/--input``: input directory for raw reads.
- ``-f/--file``: sample info table, with read pairs provided either as
  ``rawreads1``/``rawreads2`` or as an ENA/SRA ``accession``.
- ``-o/--output``: output directory.
- ``-r/--reference``: local path, URL, or NCBI genome assembly accession such
  as ``GCF_000001405.40`` for a host reference genome file.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``--fraction``: compute microbial fraction with SingleM after preprocessing.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Cataloging
^^^^^^^^^^

Assembles reads, bins contigs into MAGs, generates bin metadata, and writes
``cataloging.tsv`` with assembly, mapping, and binning summary statistics.

.. code-block:: console

   $ drakkar cataloging -i /path/to/preprocessed -o drakkar_output -m individual

Options:

- ``-i/--input``: directory with preprocessed reads or compatible workflow input.
- ``-f/--file``: sample info table. See *Read resolution* below for how the
  workflow decides which reads to use for assembly and mapping.
- ``-o/--output``: output directory.
- ``-m/--mode``: assembly modes such as ``individual`` and ``all``.
- ``-b/--binners``: comma-separated binners to run
  (``metabat``, ``maxbin``, ``semibin``, ``comebin``; default: all).
- ``--min-completeness``: minimum completeness (%) for a bin to be kept by
  Binette (default: ``70``).
- ``--max-contamination``: maximum contamination (%) for a bin to be kept by
  Binette (default: ``10``).
- ``--min-bin-length`` (``--min-length``): minimum bin length in bp
  (default: ``200000``).
- ``--max-bin-length`` (``--max-length``): maximum bin length in bp
  (default: ``10000000``).
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Read resolution when using ``-f/--file``
"""""""""""""""""""""""""""""""""""""""""

When a sample info table is provided, cataloging resolves the reads to assemble
and map in the following priority order per sample:

1. **Explicit preprocessed columns** — if the table contains
   ``preprocessedreads1`` and ``preprocessedreads2`` columns for a sample,
   those paths are used directly. Both columns must be present together.

2. **Auto-detected preprocessed reads** — if no ``preprocessedreads1``/
   ``preprocessedreads2`` columns are supplied, DRAKKAR checks whether
   ``preprocessing/final/<sample>_1.fq.gz`` and the matching R2 file exist
   inside the output directory. If they do, those quality-filtered reads are
   used. This covers the standard case of running cataloging after
   preprocessing in the same output directory without changing the input file.

3. **Raw reads or accession** — if neither of the above is available,
   cataloging falls back to ``rawreads1``/``rawreads2`` paths or an ENA/SRA
   ``accession`` from the table.

This allows a single input file to carry assembly grouping (``assembly``) and
coverage grouping (``coverage``) metadata alongside raw read paths, while
cataloging automatically upgrades to quality-filtered reads whenever they are
available.

Profiling
^^^^^^^^^

Dereplicates MAGs and maps reads to estimate abundance, with optional microbial
fraction estimation.

.. code-block:: console

   $ drakkar profiling -b /path/to/bins -R reads.tsv -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-r/--reads_dir``: directory with reads.
- ``-R/--reads_file``: sample info table with reads, using either
  ``rawreads1``/``rawreads2`` or an ENA/SRA ``accession``.
- ``-o/--output``: output directory.
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``-f/--fraction``: compute microbial fraction with SingleM.
- ``-a/--ani``: dRep ANI threshold.
- ``-n/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep.
- ``-q/--quality``: CSV/TSV with genome, completeness, and contamination; use
  this instead of CheckM2.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Genome profiling writes three tables to ``profiling_genomes/final/``:

- ``counts.tsv``: reads mapped to each dereplicated genome per sample.
- ``bases.tsv``: genome bases covered in each sample.
- ``mags.tsv``: metadata of the dereplicated genomes listed in the two tables
  above, with one row per genome and the columns ``magid``, ``completeness``,
  ``contamination``, ``size``, ``contigs``, ``longest_contig``, ``n50``,
  ``gc``, ``cluster``, ``cluster_members``, and ``score``. Completeness and
  contamination come from CheckM2 or from the file passed with ``-q/--quality``
  and are ``NA`` when ``-n/--ignore_quality`` is used without a quality file;
  ``cluster``, ``cluster_members`` (number of input bins collapsed into the
  cluster), and ``score`` come from dRep.

Annotating
^^^^^^^^^^

Annotates dereplicated MAGs taxonomically and/or functionally.
When taxonomy annotation is enabled, DRAKKAR also writes
``annotating/bacteria.tree`` and, when archaeal MAGs are present,
``annotating/archaea.tree`` by pruning GTDB-Tk classify trees down to the
input genomes only.

.. code-block:: console

   $ drakkar annotating -b /path/to/mags -o drakkar_output --annotation-type taxonomy,function

.. code-block:: console

   $ drakkar annotating -b /path/to/mags -o drakkar_output --annotation-type genes

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-o/--output``: output directory.
- ``--annotation-type``: comma-separated annotation targets:

  - ``taxonomy``: run GTDB-Tk taxonomy.
  - ``function``: run all functional components below.
  - ``genes``: run only gene-level components
    (``kegg,cazy,pfam,virulence,amr,signalp``).
  - ``kegg``: KEGG ortholog HMM annotation.
  - ``cazy``: CAZy HMM annotation.
  - ``pfam``: PFAM HMM annotation.
  - ``virulence`` (alias: ``vfdb``): VFDB-based virulence annotation.
  - ``amr``: AMR HMM annotation.
  - ``signalp``: signal peptide prediction.
  - ``dbcan``: dbCAN/CGC annotation.
  - ``antismash``: biosynthetic cluster annotation.
  - ``defense``: DefenseFinder systems and genes.
  - ``mobile`` (alias: ``genomad``): geNomad mobile and viral regions.
  - ``network``: metabolic network reconstruction.

Structure annotation with Foldseek/ProstT5 remains work in progress. It is not
an available ``--annotation-type`` in Drakkar 2.0.0, and the experimental
database preparation command is not exposed by the CLI.

- ``--gtdb-version``: GTDB release number for taxonomy annotation. DRAKKAR uses
  ``GTDB_DB_<version>`` from ``config.yaml``; if omitted, it uses ``GTDB_DB``.
- ``--annotation-evalue``: maximum fallback e-value for applicable merged gene
  annotation hits (default: ``1e-10``); sources with native cutoffs retain
  those policies.
- ``--annotation-identity``: minimum percent identity for merged gene
  annotation hits with identity values, currently VFDB/MMseqs hits
  (default: ``50``).
- ``--annotation-query-coverage``: minimum VFDB/MMseqs query coverage as a
  fraction from 0 to 1 (default: ``0.5``).
- ``--annotation-target-coverage``: minimum VFDB/MMseqs target coverage as a
  fraction from 0 to 1 (default: ``0.5``).
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Output behavior for partial functional runs:

- ``annotating/gene_annotations.tsv.xz`` is generated when any gene-level
  source is selected
  (``kegg,cazy,pfam,virulence,amr,signalp,defense``).
- ``annotating/cluster_annotations.tsv.xz`` is generated when any cluster-level
  source is selected (``dbcan,antismash,defense,mobile``).
- Merged tables are still generated from the available sources when only a
  subset of functional components is selected. Only explicitly enabled
  sources are read; leftover source files from an earlier partial run are
  ignored.
- Functional runs also write ``annotating/annotation_manifest.yaml`` and
  ``annotating/annotation_qc.tsv``.

``gene_annotations.tsv.xz`` is a long-form table with one row per retained
gene/source hit instead of one row per gene. A ``source=prodigal`` row keeps
every predicted gene in the table, while functional sources can contribute
several ranked rows. Use ``(mag, gene)`` as the gene key and
``(mag, gene, source, hit_rank)`` as the row key. ``is_primary`` marks rank 1
without discarding lower-ranked accepted evidence.

See :doc:`annotation_tables` for a before/after example, the exact ordered
schema and data types, source/method/evidence values, source-specific
``details`` fields, a complete 1.x-to-2.0 migration map, and pandas/R examples.

KOfam acceptance always requires the configured release's ``ko_list``. If the
cutoff table is missing or empty, Drakkar stops with a database-reinstallation
message rather than substituting a different scientific filter.

``cluster_annotations.tsv.xz`` contains one row per retained region or system.
Its stable columns are ``mag``, ``cluster_id``, ``contig``, ``start``, ``end``,
``source``, ``method``, ``evidence``, ``type``, ``gene_count``, ``substrate``,
``gene_functions``, ``pul_id``, and ``details``. Treat
``(mag, source, cluster_id)`` as the unique cluster key. ``details`` preserves
the complete native summary record, and source schema changes fail loudly
instead of producing blank columns.

The annotation sidecars provide run-level auditability:

- ``annotation_manifest.yaml`` records Drakkar's version, enabled sources,
  filters, configured tool modules and Conda dependencies, database paths and
  available installation checksums, plus checksums of the final annotation
  tables.
- ``annotation_qc.tsv`` reports per-MAG/source input records, retained records,
  rejected records when Drakkar performs the filtering, missing mappings or
  coordinates, unique genes/clusters, and whether filtering occurred upstream
  or during merging. A blank rejected count means the upstream tool emitted
  only accepted calls, so the rejected total is unavailable.

The 2.0 gene table is intentionally not append- or column-compatible with 1.x.
Regenerate the annotation outputs after upgrading; do not concatenate old and
new tables. The complete migration procedure, including the required VFDB
mapping rebuild, is documented under :ref:`migrating-gene-tables-2`.

Expressing
^^^^^^^^^^

Maps metatranscriptomic reads to annotated genes to quantify expression.

.. code-block:: console

   $ drakkar expressing -b /path/to/mags -R transcriptome.tsv -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-r/--reads_dir``: directory with transcriptome reads.
- ``-R/--reads_file``: transcriptome sample table, using either
  ``rawreads1``/``rawreads2`` or an ENA/SRA ``accession``.
- ``-o/--output``: output directory.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Dereplicating
^^^^^^^^^^^^^

Runs only the dereplication step and outputs dereplicated genomes to
``dereplicating/final``.

.. code-block:: console

   $ drakkar dereplicating -b /path/to/bins -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-o/--output``: output directory.
- ``-a/--ani``: dRep ANI threshold.
- ``-n/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep.
- ``-q/--quality``: CSV/TSV with genome, completeness, and contamination; use
  this instead of CheckM2.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.

Inspecting
^^^^^^^^^^

Runs microdiversity and mapping inspection workflows.

.. code-block:: console

   $ drakkar inspecting -b /path/to/mags -m /path/to/bams -c coverage.tsv -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-m/--mapping_dir``: directory with BAM files.
- ``-c/--cov_file``: coverage table per genome per sample.
- ``-o/--output``: output directory.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark`` / ``--memory-multiplier`` / ``--time-multiplier`` /
  ``--snakemake-*`` / ``--slurm-*``: see :ref:`snakemake-slurm-management`.
