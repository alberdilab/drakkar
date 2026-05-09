User guide
==========

This page introduces how DRAKKAR is organized, what kinds of inputs it expects,
and where to find the detailed workflow and operations documentation.

Quickstart
----------

Run the complete pipeline with a sample info table:

.. code-block:: console

   $ drakkar complete -f input_info.tsv -o drakkar_output

Run the complete pipeline using a directory of reads:

.. code-block:: console

   $ drakkar complete -i /path/to/reads -o drakkar_output

Core concepts
-------------

- **Modules**: DRAKKAR can be run end-to-end with ``drakkar complete`` or as
  independent modules such as preprocessing, cataloging, profiling,
  annotating, expressing, dereplicating, inspecting, database, logging,
  config, and transfer.
- **Output directory**: all outputs are written under ``-o/--output`` and
  organized into predictable module-specific folders.
- **Profiles**: use ``-p/--profile`` to select a Snakemake profile. The
  default is ``slurm``.
- **Environments**: use ``-e/--env_path`` to select a shared Conda environment
  directory.
- **Run logs**: every workflow run writes a metadata file
  ``drakkar_YYYYMMDD-HHMMSS.yaml`` and captures Snakemake stdout/stderr in
  ``log/drakkar_<run_id>.snakemake.log``.
- **Locked runs**: output-writing workflows support ``--overwrite`` to delete a
  locked output directory and rerun after a broken Snakemake session.

Input formats
-------------

You can provide inputs as read directories or as a sample info table.

Directory input
^^^^^^^^^^^^^^^

Provide a directory with paired-end reads. DRAKKAR expects matching read-pair
names such as ``*_1.fq.gz`` and ``*_2.fq.gz``.

.. code-block:: console

   $ drakkar preprocessing -i /path/to/reads -o drakkar_output

Sample info table (TSV)
^^^^^^^^^^^^^^^^^^^^^^^

A tab-separated table can include any of these columns. Only the columns needed
for the chosen workflow are required.

- ``sample``: sample name.
- ``rawreads1``: path or URL to R1 reads (raw, before preprocessing).
- ``rawreads2``: path or URL to R2 reads (raw, before preprocessing).
- ``accession``: ENA/SRA paired-end run accession such as ``ERR4303216`` or
  ``SRR12345678``. Use this instead of ``rawreads1`` and ``rawreads2`` when
  you want DRAKKAR to download the read pair automatically.
- ``preprocessedreads1``: explicit path to quality-filtered R1 reads for use
  in cataloging. Takes priority over all other read columns. See
  *Cataloging read resolution* below.
- ``preprocessedreads2``: explicit path to quality-filtered R2 reads for use
  in cataloging. Must be provided together with ``preprocessedreads1``.
- ``reference_name``: host reference label for host-removal workflows.
- ``reference_path``: local path or URL to a host FASTA, or to a tarball
  containing the FASTA plus Bowtie2 index files.
- ``assembly``: labels defining assembly groups. Legacy ``coassembly`` is
  still accepted.
- ``coverage``: labels defining coverage-sharing groups for multicoverage
  cataloging.

Example:

.. code-block:: text

   sample\trawreads1\trawreads2\taccession\treference_name\treference_path\tassembly\tcoverage
   sample1\tpath/sample1_1.fq.gz\tpath/sample1_2.fq.gz\t\tref1\tpath/ref1.fna\tassembly1,all\tcoverage1
   sample2\t\t\tERR4303216\tref1\tpath/ref1.fna\tassembly2,all\tcoverage2

Input notes
^^^^^^^^^^^

- Read files can be local paths or remote URLs (http/https/ftp).
- Sample tables can also use an ``accession`` column with ENA/SRA paired-end
  run accessions; DRAKKAR resolves and downloads the matching R1 and R2 FASTQ
  files automatically.
- ``-r/--reference``, ``-x/--reference-index``, and ``reference_path`` values
  can be local files or remote URLs.
- Reference inputs may be FASTA files, compressed FASTA files, or tarballs
  containing a FASTA plus Bowtie2 index files.
- Genome lists passed through options such as ``-B/--bins_file`` can also use
  remote URLs; DRAKKAR caches them locally before execution.
- Directory-style inputs such as ``-i/--input`` and ``-b/--bins_dir`` must be
  local filesystem paths.
- The preferred sample-table column name is ``assembly``. The legacy column
  name ``coassembly`` is still accepted.
- Assembly labels can be any identifiers you choose; they do not need to match
  sample names.
- ``-m individual`` adds per-sample assemblies alongside grouped assemblies.
- ``-b/--binners`` selects the binners used in cataloging. Use a comma-separated
  list of ``metabat``, ``maxbin``, ``semibin``, and ``comebin``; the default is
  all four.
- ``--multicoverage`` maps samples sharing the same coverage label to each
  other's individual assemblies.

Cataloging read resolution
^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``drakkar cataloging`` (or ``drakkar complete``) loads a sample info table
with ``-f/--file``, it resolves the reads to use for assembly and mapping in the
following priority order for each sample:

1. **``preprocessedreads1`` / ``preprocessedreads2`` columns** — if both are
   present the cataloging workflow uses these paths directly. This is the
   explicit override for cases where preprocessed reads live outside the
   default output tree.

2. **``preprocessing/final/<sample>_1.fq.gz``** — if neither
   ``preprocessedreads1`` nor ``preprocessedreads2`` is supplied but a prior
   ``drakkar preprocessing`` run has already written quality-filtered reads into
   the output directory, cataloging detects and uses them automatically. This
   is the typical case when running cataloging as a follow-up step after
   preprocessing in the same output directory.

3. **``rawreads1`` / ``rawreads2`` or ``accession``** — fallback to raw input
   paths. This path is taken when neither preprocessed column is present and no
   ``preprocessing/final/`` files are found. The assembly will run directly on
   unfiltered reads.

This means you can keep a single input table that contains raw read paths (or
accessions) together with ``assembly`` and ``coverage`` grouping columns, and
cataloging will automatically pick up the quality-filtered reads from a
completed preprocessing run without any changes to the file.

Guide map
---------

Use the next pages depending on what you need:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Topic
     - Where to go next
   * - Running the complete workflow or a specific module
     - See :doc:`workflows`.
   * - Databases, logging, config, transfer, outputs, and troubleshooting
     - See :doc:`operations`.
   * - Command list only
     - See :doc:`api`.
