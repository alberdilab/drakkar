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
  cataloging, profiling, annotating, expressing, dereplicating, database,
  config, inspecting, transfer).
- **Output directory**: All outputs are written under ``-o/--output``.
- **Profiles**: Use ``-p/--profile`` to select a Snakemake profile (default:
  ``slurm``).
- **Environments**: Use ``-e/--env_path`` to set a shared Conda env directory.
- **Run logs**: Every run writes a metadata file ``drakkar_YYYYMMDD-HHMMSS.yaml``
  to the output directory.
- **Locked runs**: Output-writing workflows support ``--overwrite`` to delete a
  locked output directory and rerun from scratch after a broken Snakemake
  session.

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
- ``reference_path``: local path or URL to a host reference FASTA, or to a
  tarball containing the FASTA plus Bowtie2 index files.
- ``coassembly``: labels that define co-assembly groups.
- ``coverage``: mapping coverage groups for multicoverage.

Example:

.. code-block:: text

   sample\trawreads1\trawreads2\treference_name\treference_path\tcoassembly\tcoverage
   sample1\tpath/sample1_1.fq.gz\tpath/sample1_2.fq.gz\tref1\tpath/ref1.fna\tassembly1,all\tcoverage1
   sample2\tpath/sample2_1.fq.gz\tpath/sample2_2.fq.gz\tref1\tpath/ref1.fna\tassembly2,all\tcoverage2

Notes:

- Read files can be local paths or remote URLs (http/https/ftp).
- ``-r/--reference``, ``-x/--reference-index``, and ``reference_path`` values
  can point to either local files or remote URLs (http/https/ftp). Reference
  inputs may be FASTA/compressed FASTA files or tarballs containing the FASTA
  plus Bowtie2 index files. Genome path lists provided through options such as
  ``-B/--bins_file`` can also use remote URLs; DRAKKAR downloads them into a
  local cache before running.
- Directory-style inputs such as ``-i/--input`` and ``-b/--bins_dir`` must remain
  local filesystem paths.
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
- ``-r/--reference``: local path or URL to a host reference genome for preprocessing.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``-m/--mode``: assembly modes (``individual``, ``all``).
- ``-t/--type``: profiling type (``genomes`` or ``pangenomes``).
- ``--annotation-type``: comma-separated targets. See *Annotating* for all values.
- ``-c/--multicoverage``: enable multicoverage mapping.
- ``--fraction``: compute microbial fraction with singlem.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil during preprocessing.
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
- ``-r/--reference``: local path or URL to a host reference genome file.
- ``-x/--reference-index``: local path or URL to a tarball containing a host
  reference FASTA and Bowtie2 index files; incompatible with ``-r/--reference``.
- ``--fraction``: compute microbial fraction with SingleM after preprocessing.
- ``--nonpareil``: estimate metagenomic coverage and diversity with Nonpareil after preprocessing.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Cataloging
^^^^^^^^^^

Assembles reads, bins contigs into MAGs, generates bin metadata, and writes
``cataloging.tsv`` with assembly, mapping, and binning summary statistics.

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
- ``-n/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep (skip CheckM quality filtering).
- ``-q/--quality``: CSV/TSV with genome, completeness, contamination; uses provided quality file instead of CheckM2.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Annotating
^^^^^^^^^^

Annotates dereplicated MAGs taxonomically and/or functionally.

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
  - ``genes``: run only gene-level components (``kegg,cazy,pfam,virulence,amr,signalp``).
  - ``kegg``: KEGG ortholog HMM annotation.
  - ``cazy``: CAZy HMM annotation.
  - ``pfam``: PFAM HMM annotation.
  - ``virulence`` (alias: ``vfdb``): VFDB-based virulence annotation.
  - ``amr``: AMR HMM annotation.
  - ``signalp``: signal peptide prediction.
  - ``dbcan``: dbCAN/CGC annotation (cluster-level).
  - ``antismash``: BGC cluster annotation.
  - ``defense``: DefenseFinder systems and genes.
  - ``mobile`` (alias: ``genomad``): geNomad viral/mobile regions.
  - ``network``: metabolic network reconstruction.
- ``--gtdb-version``: GTDB release number for taxonomy annotation. DRAKKAR
  uses ``GTDB_DB_<version>`` from ``config.yaml``; if omitted, it uses
  ``GTDB_DB``.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Output behavior for partial functional runs:

- ``annotating/gene_annotations.tsv.xz`` is generated when any gene-level source is selected
  (``kegg,cazy,pfam,virulence,amr,signalp,defense``).
- ``annotating/cluster_annotations.tsv.xz`` is generated when any cluster-level source is selected
  (``dbcan,antismash,defense,mobile``).
- You can run any subset of functional components; merged tables are still generated from the available sources.

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
- ``-n/--ignore_quality``: pass ``--ignoreGenomeQuality`` to dRep (skip CheckM quality filtering).
- ``-q/--quality``: CSV/TSV with genome, completeness, contamination; uses provided quality file instead of CheckM2.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Database
^^^^^^^^

Installs or updates one managed annotation database release at a time. This is
a maintenance workflow and is not triggered by ``drakkar complete``.

Supported database subcommands:

- ``kegg`` (alias: ``kofams``)
- ``cazy``
- ``pfam``
- ``vfdb``
- ``amr``

.. code-block:: console

   $ drakkar database amr --directory /projects/alberdilab/data/databases/drakkar/amr --version 2025-07-16.1

.. code-block:: console

   $ drakkar database kegg --directory /projects/alberdilab/data/databases/drakkar/kofams --version 2026-02-01 --set-default

.. code-block:: console

   $ drakkar database kegg --directory /projects/alberdilab/data/databases/drakkar/kofams --version 2026-02-01 --download-runtime 180

.. code-block:: console

   $ drakkar database cazy --directory /projects/alberdilab/data/databases/drakkar/cazy --version V14 --set-default

.. code-block:: console

   $ drakkar database pfam --directory /projects/alberdilab/data/databases/drakkar/pfam --version Pfam37.4 --set-default

.. code-block:: console

   $ drakkar database vfdb --directory /projects/alberdilab/data/databases/drakkar/vfdb --set-default

Options:

- ``--directory``: base directory where the release folder will be created.
- ``--version``: folder name to create inside ``--directory``. For ``kegg``,
  use the KEGG archive date such as ``2026-02-01``. For ``cazy``, use the
  upstream dbCAN release label such as ``V14``. For ``pfam``, use the Pfam
  release directory name such as ``Pfam37.4``. For ``amr``, use the NCBI
  AMRFinder release directory name such as ``2025-07-16.1``. For ``vfdb``,
  you can omit ``--version`` and DRAKKAR will use the download date in UTC,
  for example ``2026-04-24``.
- ``--download-runtime``: runtime in minutes for the database
  download/preparation rule (default: ``120``).
- ``--set-default``: update the corresponding database path in ``config.yaml``
  after the installation finishes successfully.
- ``-e/--env_path``: shared Conda env dir.
- ``-p/--profile``: Snakemake profile.

Behavior:

- The selected database is installed into ``--directory/--version/``.
- For managed annotation databases, ``config.yaml`` stores the release
  directory, not the internal HMM or MMseqs prefix file.
- The workflow resolves the expected internal files automatically, for example
  ``kofams``, ``pfam``, ``amr.tsv``, or ``vfdb``.
- ``--set-default`` rewrites that config entry to the newly installed release
  directory.

Database-specific rules:

- ``kegg`` (alias: ``kofams``): use a KEGG archive date in ``YYYY-MM-DD``
  format, such as ``2026-02-01``. DRAKKAR downloads ``profiles.tar.gz`` from
  ``https://www.genome.jp/ftp/db/kofam/archives/<version>/``, extracts the HMM
  profiles, concatenates them into a single ``kofams`` database, downloads the
  KEGG hierarchy JSON, and runs ``hmmpress``. If the archive is missing,
  DRAKKAR points you to ``https://www.genome.jp/ftp/db/kofam/archives/``. The
  default ``--download-runtime`` is ``120`` minutes and is mainly intended for
  this large download.
- ``cazy``: use the dbCAN release label, such as ``V14``. DRAKKAR downloads the
  dbCAN HMM database from
  ``https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/<version>/dbCAN-HMMdb-<version>.txt``
  and runs ``hmmpress``. If the requested release is missing, DRAKKAR points
  you to ``https://pro.unl.edu/dbCAN2/browse_download.php``.
- ``pfam``: use the Pfam release directory name, such as ``Pfam37.4``. DRAKKAR
  downloads ``Pfam-A.hmm.gz`` from
  ``https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/<version>/``, downloads
  the EC mapping table, unzips the HMM file, and runs ``hmmpress``. If the
  requested release is missing, DRAKKAR points you to
  ``https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/``.
- ``amr``: use the NCBI AMRFinder release directory name, such as
  ``2025-07-16.1``. DRAKKAR downloads both
  ``NCBIfam-AMRFinder.HMM.tar.gz`` and ``NCBIfam-AMRFinder.tsv`` from
  ``https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/<version>/``, merges the
  extracted HMMs into one database, and runs ``hmmpress``. If the requested
  release is missing, DRAKKAR points you to
  ``https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/``.
- ``vfdb``: there is no upstream version directory. DRAKKAR downloads the
  current ``VFDB_setB_pro.fas.gz`` from
  ``https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz``, creates the MMseqs2
  database, and if ``--version`` is omitted it uses the UTC download date as
  the release folder and logged version, for example ``2026-04-24``.

Version logging:

- Each run writes ``database_versions.yaml`` inside the installed release directory.
- The log records the requested version, resolved install directory, source URLs,
  source-version label, and SHA256 checksums and file sizes for the installed
  assets.

Config
^^^^^^

Views or edits the installed DRAKKAR configuration file at
``drakkar/workflow/config.yaml``.

.. code-block:: console

   $ drakkar config --view

.. code-block:: console

   $ drakkar config --edit

Options:

- ``--view``: print the config file path and contents.
- ``--edit``: open the config file in a terminal editor.

Behavior:

- ``--edit`` uses ``$VISUAL``, then ``$EDITOR``, then falls back to ``nano``,
  ``vim``, or ``vi``.
- The command edits the installed package config directly, so changes affect
  later workflow runs from that installation.

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
- ``cataloging.tsv`` assembly, mapping, and binning summary table.
- ``profiling_genomes/`` dereplication, mapping, and abundance tables.
- ``profiling_pangenomes/`` pangenome profiling outputs.
- ``annotating/`` annotation tables.
- ``expressing/`` expression outputs.
- ``dereplicating/`` dereplicated genomes (dereplicating-only mode).
- ``<directory>/<version>/database_versions.yaml`` installation log for a managed database release.

Troubleshooting
---------------

- **Locked directory**: run ``drakkar unlock -o <output_dir>``.
- **Missing bins**: provide ``-b/--bins_dir`` or ``-B/--bins_file``.
- **Missing reads**: provide ``-r/--reads_dir`` or ``-R/--reads_file``.
- **SFTP errors**: ensure the remote directory exists and credentials are valid.
