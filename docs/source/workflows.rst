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
- ``-f/--file``: sample info table (TSV), with read pairs provided either as
  ``rawreads1``/``rawreads2`` or as an ENA/SRA ``accession``.
- ``-o/--output``: output directory.
- ``-r/--reference``: local path or URL to a host reference genome for preprocessing.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``-m/--mode``: assembly modes such as ``individual`` and ``all``.
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``--annotation-type``: comma-separated annotation targets. See *Annotating*
  below for the full set.
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``--fraction``: compute microbial fraction with SingleM.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil.
- ``-a/--ani``: dRep ANI threshold (default: ``0.98``).
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile (default: ``slurm``).

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
- ``-r/--reference``: local path or URL to a host reference genome file.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``--fraction``: compute microbial fraction with SingleM after preprocessing.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.

Cataloging
^^^^^^^^^^

Assembles reads, bins contigs into MAGs, generates bin metadata, and writes
``cataloging.tsv`` with assembly, mapping, and binning summary statistics.

.. code-block:: console

   $ drakkar cataloging -i /path/to/preprocessed -o drakkar_output -m individual

Options:

- ``-i/--input``: directory with preprocessed reads or compatible workflow input.
- ``-f/--file``: sample info table, with read pairs provided either as
  ``rawreads1``/``rawreads2`` or as an ENA/SRA ``accession``.
- ``-o/--output``: output directory.
- ``-m/--mode``: assembly modes such as ``individual`` and ``all``.
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.

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
- ``--gtdb-version``: GTDB release number for taxonomy annotation. DRAKKAR uses
  ``GTDB_DB_<version>`` from ``config.yaml``; if omitted, it uses ``GTDB_DB``.
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile.

Output behavior for partial functional runs:

- ``annotating/gene_annotations.tsv.xz`` is generated when any gene-level
  source is selected
  (``kegg,cazy,pfam,virulence,amr,signalp,defense``).
- ``annotating/cluster_annotations.tsv.xz`` is generated when any cluster-level
  source is selected (``dbcan,antismash,defense,mobile``).
- Merged tables are still generated from the available sources when only a
  subset of functional components is selected.

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
