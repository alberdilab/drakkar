Usage
=====

This page describes how to run DRAKKAR, the expected inputs, and the outputs
produced by each module.

Quickstart
----------

Run the complete pipeline with a sample info table:

.. code-block:: console

   $ drakkar complete -f input_info.tsv -o drakkar_output

Run the complete pipeline using a directory of reads:

.. code-block:: console

   $ drakkar complete -i /path/to/reads -o drakkar_output

Key concepts
------------

- **Modules**: Each workflow stage can be run independently (preprocessing,
  cataloging, profiling, annotating, expressing, dereplicating, inspecting).
- **Output directory**: All outputs are written under ``-o/--output``.
- **Profiles**: Use ``-p/--profile`` to select a Snakemake profile (default:
  ``slurm``).
- **Environments**: Use ``-e/--env_path`` to set a shared Conda env directory.
- **Run logs**: Every run writes a metadata file ``drakkar_YYYYMMDD-HHMMSS.yaml``
  to the output directory.

Input formats
-------------

You can provide inputs as a directory of reads or as a sample info table.

Directory input
^^^^^^^^^^^^^^^

Provide a directory with paired-end reads. DRAKKAR expects read pairs with
consistent naming conventions (e.g., ``*_1.fq.gz`` and ``*_2.fq.gz``).

.. code-block:: console

   $ drakkar preprocessing -i /path/to/reads -o drakkar_output

Sample info table (TSV)
^^^^^^^^^^^^^^^^^^^^^^^

A tab-separated table can include any of these columns (only those required
by the selected module are needed):

- ``sample``: sample name (required).
- ``rawreads1``: path or URL to R1 reads.
- ``rawreads2``: path or URL to R2 reads.
- ``reference_name``: host reference label (preprocessing with host removal).
- ``reference_path``: host reference fasta path.
- ``coassembly``: labels that define co-assembly groups.
- ``coverage``: mapping coverage groups for multicoverage.

Example:

.. code-block:: text

   sample\trawreads1\trawreads2\treference_name\treference_path\tcoassembly\tcoverage
   sample1\tpath/sample1_1.fq.gz\tpath/sample1_2.fq.gz\tref1\tpath/ref1.fna\tassembly1,all\tcoverage1
   sample2\tpath/sample2_1.fq.gz\tpath/sample2_2.fq.gz\tref1\tpath/ref1.fna\tassembly2,all\tcoverage2

Notes:

- Read files can be local paths or remote URLs (http/https/ftp).
- The ``coassembly`` column groups samples for pooled assemblies.
- ``-m individual`` adds per-sample assemblies in addition to co-assemblies.
- ``--multicoverage`` maps samples that share the same coverage label to each
  other's individual assemblies (not compatible with co-assemblies).

Complete workflow
-----------------

Run the full pipeline in sequence:

.. code-block:: console

   $ drakkar complete -f input_info.tsv -o drakkar_output -m individual -t genomes

Options:

- ``-i/--input``: input directory for reads.
- ``-f/--file``: sample info table (TSV).
- ``-o/--output``: output directory.
- ``-r/--reference``: host reference genome for preprocessing.
- ``-m/--mode``: assembly modes (``individual``, ``all``).
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``--annotation-type``: ``taxonomy``, ``function``, or both.
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``--fraction``: compute microbial fraction with singlem.
- ``-a/--ani``: dRep ANI threshold (default: 0.98).
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile (default: slurm).

Module reference
----------------

Preprocessing
^^^^^^^^^^^^^

Quality filters reads, optionally removes host DNA, and writes cleaned reads
and summary tables.

.. code-block:: console

   $ drakkar preprocessing -i /path/to/reads -o drakkar_output -r host.fna

Options:

- ``-i/--input``: input directory for raw reads.
- ``-f/--file``: sample info table.
- ``-o/--output``: output directory.
- ``-r/--reference``: host reference genome file.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Cataloging
^^^^^^^^^^

Assembles reads, bins contigs into MAGs, and generates bin metadata.

.. code-block:: console

   $ drakkar cataloging -i /path/to/preprocessed -o drakkar_output -m individual

Options:

- ``-i/--input``: directory with preprocessed reads or sample info output.
- ``-f/--file``: sample info table.
- ``-o/--output``: output directory.
- ``-m/--mode``: assembly modes (``individual``, ``all``).
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``-e/--env_path``: shared Conda env dir.
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
- ``-R/--reads_file``: sample info table with reads.
- ``-o/--output``: output directory.
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``-f/--fraction``: compute microbial fraction with singlem.
- ``-a/--ani``: dRep ANI threshold.
- ``-q/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep (skip CheckM quality filtering).
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Annotating
^^^^^^^^^^

Annotates dereplicated MAGs taxonomically and/or functionally.

.. code-block:: console

   $ drakkar annotating -b /path/to/mags -o drakkar_output --annotation-type taxonomy,function

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-o/--output``: output directory.
- ``--annotation-type``: ``taxonomy``, ``function``, or both.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Expressing
^^^^^^^^^^

Maps metatranscriptomic reads to annotated genes to quantify expression.

.. code-block:: console

   $ drakkar expressing -b /path/to/mags -R transcriptome.tsv -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-r/--reads_dir``: directory with transcriptome reads.
- ``-R/--reads_file``: transcriptome sample table.
- ``-o/--output``: output directory.
- ``-e/--env_path``: shared Conda env dir.
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
- ``-q/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep (skip CheckM quality filtering).
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Inspecting
^^^^^^^^^^

Runs microdiversity and mapping inspection workflows.

.. code-block:: console

   $ drakkar inspecting -b /path/to/mags -m /path/to/bams -c coverage.tsv -o drakkar_output

Options:

- ``-b/--bins_dir``: directory with MAG/bin FASTA files.
- ``-B/--bins_file``: file listing MAG/bin paths.
- ``-m/--mapping_dir``: directory with mapping (BAM) files.
- ``-c/--cov_file``: coverage table per genome per sample.
- ``-o/--output``: output directory.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Transfer
^^^^^^^^

Transfers outputs via SFTP while preserving the original folder structure.
The remote base directory must exist.

.. code-block:: console

   $ drakkar transfer --host example.org --user you -l drakkar_output -r /remote/path --results -v

Flags:

- ``--all``: transfer the entire output directory.
- ``--data``: transfer everything except ``.snakemake``.
- ``--results``: transfer the union of ``-a/-m/-p/-b/-e``.
- ``-a/--annotations``: annotation outputs.
- ``-m/--mags``: dereplicated MAGs.
- ``-p/--profile``: profiling outputs.
- ``-e/--expression``: expression outputs.
- ``-b/--bins``: cataloging bins (recursive).
- ``--erda``: use ERDA defaults (``io.erda.dk``).
- ``-v/--verbose``: log each transfer.

Maintenance commands
--------------------

Unlock a working directory if Snakemake left a lock:

.. code-block:: console

   $ drakkar unlock -o drakkar_output

Update Drakkar in the current environment:

.. code-block:: console

   $ drakkar update

Outputs
-------

Key output locations:

- ``preprocessing/`` cleaned reads and summaries.
- ``cataloging/`` assemblies, bins, and metadata.
- ``profiling_genomes/`` dereplication, mapping, and abundance tables.
- ``profiling_pangenomes/`` pangenome profiling outputs.
- ``annotating/`` annotation tables.
- ``expressing/`` expression outputs.
- ``dereplicating/`` dereplicated genomes (dereplicating-only mode).

Troubleshooting
---------------

- **Locked directory**: run ``drakkar unlock -o <output_dir>``.
- **Missing bins**: provide ``-b/--bins_dir`` or ``-B/--bins_file``.
- **Missing reads**: provide ``-r/--reads_dir`` or ``-R/--reads_file``.
- **SFTP errors**: ensure the remote directory exists and credentials are valid.
