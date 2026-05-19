Operations and troubleshooting
==============================

This page covers the operational commands in DRAKKAR: database preparation,
configuration, status inspection, logging, result transfer, output layout, and
common recovery tasks.

Operations overview
-------------------

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Command
     - Purpose
     - Typical use
   * - ``drakkar database``
     - Install or update supported annotation database releases.
     - Prepare KEGG, CAZy, PFAM, AMR, or VFDB resources before annotation.
   * - ``drakkar config``
     - View or edit the installed workflow configuration.
     - Inspect or change database paths and default settings.
   * - ``drakkar logging``
     - Inspect workflow metadata and Snakemake logs.
     - Diagnose failed runs, locked directories, and progress state.
   * - ``drakkar status``
     - Show rule and sample progress for a workflow run.
     - Monitor the latest run, a selected output directory, or one metadata YAML.
   * - ``drakkar transfer``
     - Transfer selected outputs by SFTP while preserving structure.
     - Move results from cluster storage to long-term or collaborator storage.
   * - ``drakkar unlock``
     - Remove a Snakemake lock from a broken output directory.
     - Recover after interrupted runs.
   * - ``drakkar update``
     - Reinstall DRAKKAR from the Git repository in the current environment.
     - Refresh the installed CLI and workflow package.

See also :ref:`snakemake-slurm-management` for the Snakemake and SLURM override
flags available on every workflow command.

Database
--------

Installs or updates one managed annotation database release at a time. This is
a maintenance workflow and is not triggered by ``drakkar complete``.

Supported database subcommands:

- ``kegg`` (alias: ``kofams``)
- ``cazy``
- ``pfam``
- ``vfdb``
- ``amr``

Examples:

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
  you can omit ``--version`` and DRAKKAR will use the UTC download date.
- ``--download-runtime``: runtime in minutes for the database download and
  preparation rule (default: ``120``).
- ``--set-default``: update the corresponding database path in ``config.yaml``
  after installation.
- ``-e/--env_path``: shared Conda environment directory.
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
  the release folder and logged version.

Version logging:

- Each run writes ``database_versions.yaml`` inside the installed release
  directory.
- The log records the requested version, resolved install directory, source
  URLs, source-version label, and installed asset checksums and file sizes.

Config
------

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

.. _snakemake-slurm-management:

Snakemake and SLURM management
-------------------------------

Every workflow subcommand (``complete``, ``preprocessing``, ``cataloging``,
``profiling``, ``annotating``, ``expressing``, ``dereplicating``,
``inspecting``, ``database``, and ``environments``) accepts the flags
described in this section. They let you tune resource limits, override
Snakemake profile settings, and pass SLURM directives without editing profile
files.

Resource caps (config.yaml)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``drakkar/workflow/config.yaml`` contains four resource-related keys that act
as cluster-wide guardrails:

- ``SNAKEMAKE_MAX_GB``: maximum memory any single rule may request, in
  gigabytes. Default: ``1024``. Dynamic per-rule memory requests are capped at
  this value.
- ``SNAKEMAKE_MAX_TIME``: maximum runtime any single rule may request, in
  minutes. Default: ``20160`` (14 days).
- ``MEMORY_MULTIPLIER``: a global integer factor applied to every per-rule
  memory request before the ``SNAKEMAKE_MAX_GB`` cap is enforced. Default:
  ``1``. Increase this when a workflow consistently runs out of memory due to
  unusually large samples.
- ``TIME_MULTIPLIER``: equivalent factor for runtime requests before the
  ``SNAKEMAKE_MAX_TIME`` cap. Default: ``1``. Increase when jobs time out on a
  slow or heavily loaded cluster.

Edit these values with ``drakkar config --edit`` or set them on the command
line with the flags below.

Resource multiplier flags
^^^^^^^^^^^^^^^^^^^^^^^^^

``--memory-multiplier N`` and ``--time-multiplier N`` apply the same scaling
as ``MEMORY_MULTIPLIER`` / ``TIME_MULTIPLIER`` in ``config.yaml`` but without
permanently changing the installed config. The command-line value overrides the
config value for that run only.

.. code-block:: console

   $ drakkar cataloging -f input.tsv -o drakkar_output --memory-multiplier 2

   $ drakkar profiling -b /path/to/bins -o drakkar_output --time-multiplier 3

Both flags accept any positive integer. They are most useful when a specific
workflow run is expected to be unusually resource-intensive.

Snakemake override flags
^^^^^^^^^^^^^^^^^^^^^^^^

These flags override the corresponding settings in the active Snakemake
profile without modifying profile files. All are optional; omitting a flag
leaves the profile value in effect.

- ``--snakemake-jobs N``: maximum number of concurrent SLURM jobs. Overrides
  the profile value (typical default: ``100``).
- ``--snakemake-cores N``: maximum local CPU cores when using the local
  executor. Overrides the profile value.
- ``--snakemake-executor EXECUTOR``: Snakemake executor plugin, e.g. ``slurm``
  or ``local``. Overrides the profile value.
- ``--snakemake-latency-wait N``: seconds to wait for output files before
  failing a rule. Overrides the profile value (slurm default: ``300``, local
  default: ``60``). Raise this on shared filesystems with high metadata
  latency.
- ``--snakemake-retries N``: number of times to retry a failed job. Overrides
  the profile value (slurm default: ``3``).
- ``--snakemake-rerun-incomplete``: force rerun of jobs whose output files were
  left incomplete by a previous interrupted run.
- ``--snakemake-keep-going``: continue running independent jobs after a failure
  instead of stopping immediately.

Examples:

.. code-block:: console

   $ drakkar complete -f input.tsv -o drakkar_output --snakemake-jobs 50 --snakemake-retries 5

   $ drakkar cataloging -f input.tsv -o drakkar_output --snakemake-executor local --snakemake-cores 32

   $ drakkar profiling -b bins/ -o drakkar_output --snakemake-rerun-incomplete --snakemake-keep-going

SLURM override flags
^^^^^^^^^^^^^^^^^^^^

These flags inject SLURM directives into Snakemake's ``--default-resources``
without requiring changes to the SLURM profile or cluster config.

- ``--slurm-partition NAME``: SLURM partition (queue) to submit all jobs to.
- ``--slurm-account NAME``: SLURM billing account.
- ``--slurm-constraint EXPR``: node constraint expression, e.g. ``gpu`` or
  ``skylake``.
- ``--slurm-nodes N``: number of nodes per SLURM job (default: ``1``).
- ``--slurm-nodelist NODES``: restrict jobs to a specific node or node list,
  e.g. ``node01`` or ``node[01-03]``.
- ``--slurm-extra ARGS``: arbitrary extra ``sbatch`` arguments passed verbatim,
  e.g. ``'--mail-type=END --mail-user=you@example.com'``.

Examples:

.. code-block:: console

   $ drakkar complete -f input.tsv -o drakkar_output --slurm-partition gpu --slurm-account myproject

   $ drakkar annotating -b bins/ -o drakkar_output --slurm-extra '--mail-type=END --mail-user=you@example.com'

SLURM benchmarking
^^^^^^^^^^^^^^^^^^

After each workflow run, DRAKKAR queries ``sacct`` for the jobs submitted
during that run and writes a resource-efficiency summary. This produces:

- ``benchmark/``: per-job resource tables under the output directory.
- ``drakkar_<run_id>_resources.yaml``: root-level summary of CPU time, memory
  peaks, and efficiency ratios for the run.

The resource summary is also shown by ``drakkar logging`` alongside the
workflow execution summary.

To skip benchmark collection, pass ``--skip-benchmark`` to any workflow
command:

.. code-block:: console

   $ drakkar preprocessing -i /path/to/reads -o drakkar_output --skip-benchmark

Status
------

Shows progress for the latest or selected Drakkar workflow run without
restarting Snakemake.

.. code-block:: console

   $ drakkar status

.. code-block:: console

   $ drakkar status -d drakkar_output --rules

.. code-block:: console

   $ drakkar status drakkar_20260510-032711.yaml --samples

Options:

- ``target``: optional output directory or ``drakkar_<run_id>.yaml`` metadata
  file. If omitted, DRAKKAR inspects the current directory.
- ``-d/--directory`` or ``-o/--output``: output directory to inspect.
- ``--run``: specific run ID or ``drakkar_<run_id>.yaml`` file name.
- ``--rules``: show rule-focused progress only.
- ``--samples``: show sample-focused progress only.
- ``--complete``: include helper rules that are hidden by default.

Behavior:

- The default view shows overall progress, rule progress for main rules, and
  sample-stage progress.
- Rule totals are parsed from the captured Snakemake job stats and completion
  lines in ``log/drakkar_<run_id>.snakemake.log``.
- Sample stages are inferred from observed sample or assembly wildcards and the
  workflow sample dictionaries under ``data/``.

Logging
-------

Inspects workflow metadata and persistent Snakemake logs to troubleshoot failed
or interrupted runs.

.. code-block:: console

   $ drakkar logging -o drakkar_output

.. code-block:: console

   $ drakkar logging -o drakkar_output --summary

.. code-block:: console

   $ drakkar logging -o drakkar_output --run 20260503-101530 --paths

Options:

- ``-o/--output``: output directory to inspect.
- ``--run``: specific run ID (``YYYYMMDD-HHMMSS``) or
  ``drakkar_<run_id>.yaml`` file name.
- ``--summary``: print only the parsed workflow summary.
- ``--tail``: number of trailing log lines to show if no failure excerpt is
  found and ``--summary`` is not used (default: ``50``).
- ``--full``: print the full Snakemake log.
- ``--paths``: list relevant metadata and log file paths.
- ``--list``: list available workflow runs in the output directory.

Behavior:

- Workflow runs write root metadata files such as
  ``drakkar_20260503-101530.yaml``.
- Snakemake stdout/stderr is captured persistently in
  ``log/drakkar_20260503-101530.snakemake.log``.
- The default logging view includes a parsed execution summary with planned
  jobs, observed rule executions, workflow progress, and detected error types.
- If the output directory is locked, run ``drakkar logging -o <output_dir>``
  before using ``drakkar unlock`` or ``--overwrite``.

Transfer
--------

Transfers outputs via SFTP while preserving the original folder structure.
The remote base directory must already exist.

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
- ``-b/--bins``: cataloging bins recursively.
- ``--erda``: use ERDA defaults (``io.erda.dk``).
- ``-v/--verbose``: log each transfer.

Maintenance commands
--------------------

Unlock a working directory if Snakemake left a lock:

.. code-block:: console

   $ drakkar unlock -o drakkar_output

Update DRAKKAR in the current environment:

.. code-block:: console

   $ drakkar update

Pass ``--skip-deps`` to refresh the package without reinstalling Python
dependencies (useful when only the workflow scripts have changed):

.. code-block:: console

   $ drakkar update --skip-deps

Outputs
-------

Key output locations:

- ``preprocessing/``: cleaned reads and preprocessing summaries.
- ``cataloging/``: assemblies, bins, and bin metadata.
- ``cataloging.tsv``: assembly, mapping, and binning summary table.
- ``profiling_genomes/``: dereplication, mapping, and abundance tables.
- ``profiling_pangenomes/``: pangenome profiling outputs.
- ``annotating/``: annotation tables.
- ``expressing/``: expression outputs.
- ``dereplicating/``: dereplicated genomes in dereplication-only mode.
- ``benchmark/``: per-SLURM-job resource tables written after each workflow run.
- ``drakkar_<run_id>.yaml``: workflow run metadata.
- ``drakkar_<run_id>_resources.yaml``: root-level SLURM resource-efficiency
  summary for the run (CPU time, memory peaks, and efficiency ratios).
- ``log/drakkar_<run_id>.snakemake.log``: persistent Snakemake stdout/stderr
  capture for a workflow run.
- ``<directory>/<version>/database_versions.yaml``: installation log for a
  managed database release.

Troubleshooting
---------------

- **Locked directory**: first run ``drakkar logging -o <output_dir>`` to inspect
  the latest workflow log, then use ``drakkar unlock -o <output_dir>`` or rerun
  with ``--overwrite``.
- **Missing bins**: provide ``-b/--bins_dir`` or ``-B/--bins_file``.
- **Missing reads**: provide ``-r/--reads_dir`` or ``-R/--reads_file``.
- **SFTP errors**: ensure the remote directory exists and the credentials are
  valid.
