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
     - Install supported annotation database releases, check the configured ones
       against their sources with ``drakkar database latest``, or bring them all
       up to date with ``drakkar database update``.
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
- ``latest`` (reports newer releases instead of installing one; see
  `Checking for newer releases`_)
- ``update`` (installs every outdated release in one command; see
  `Updating every outdated database`_)

Foldseek/ProstT5 structure annotation and its database installer remain work in
progress and are not available through the Drakkar 2.0.0 CLI.

Examples:

.. code-block:: console

   $ drakkar database amr --version 2025-07-16.1

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
  Optional: it defaults to ``DATABASES_DIR/<database>`` from ``config.yaml``,
  so ``drakkar database pfam --version Pfam38.2`` installs into
  ``DATABASES_DIR/pfam/Pfam38.2``. KEGG releases live under ``kofams``. Pass
  ``--directory`` to install somewhere else, and set ``DATABASES_DIR`` with
  ``drakkar config --edit`` if the command reports that it is unset.
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
  the release folder and logged version. Drakkar 2.0 uses a corrected,
  schema-marked mapping table; install a fresh dated release after upgrading
  because 1.x mappings can label the organism as the virulence-factor type.

Version logging:

- Each run writes ``database_versions.yaml`` inside the installed release
  directory.
- The log records the requested version, resolved install directory, source
  URLs, source-version label, and installed asset checksums and file sizes.

Checking for newer releases
---------------------------

``drakkar database latest`` queries each database source and compares its newest
release with the version wired into ``config.yaml``. Nothing is downloaded or
installed, and no output directory is touched.

.. code-block:: console

   $ drakkar database latest

.. code-block:: console

   $ drakkar database latest kegg pfam

.. code-block:: text

   DATABASE               CONFIGURED       LATEST           STATUS
   amr                    2025-07-16.1     2026-08-07.1     outdated
   cazy                   V14              V15              outdated
   gtdb                   232              232              up to date

Every outdated database is followed by the command that installs the newer
release, so the reported version can be applied directly.

Databases that can be checked:

- ``kegg`` (alias: ``kofams``), ``cazy``, ``pfam``, ``vfdb``, ``amr`` and
  ``foldseek``, using the release directory recorded in ``config.yaml``.
- ``gtdb``, read from the ``GTDB_DB`` entry. GTDB reference data is installed by
  GTDB-Tk rather than by DRAKKAR, so its row reports the newest release without
  offering an install command.

Options:

- ``databases``: names to check. All of them are checked when none is given.
- ``--timeout``: seconds to wait for each source before giving up
  (default: ``20``).

Behavior:

- ``vfdb`` and ``foldseek`` have no upstream version numbers, so their releases
  are compared by date: the download date the release folder is named after
  against the ``Last-Modified`` date of the source file.
- ``cazy`` publishes no directory listing, so consecutive dbCAN release numbers
  are requested until one is missing.
- A source that cannot be reached is reported as ``unknown`` for that database
  only. The command still reports every other database and exits ``0``, so an
  unreachable mirror never blocks the check.
- Installing a newer release does not rebuild existing outputs. See
  `Cross-run consistency check`_ for what happens on the next run in a
  directory built with the earlier release.

Updating every outdated database
--------------------------------

``drakkar database update`` runs the version check above and then installs the
newest release of every managed database that is behind, one after another,
without naming each one by hand.

.. code-block:: console

   $ drakkar database update

.. code-block:: console

   $ drakkar database update --yes

.. code-block:: console

   $ drakkar database update kegg pfam --yes

The command prints its plan and stops unless ``--yes`` is given, because a full
update downloads tens of gigabytes and then repoints ``config.yaml`` at the new
releases:

.. code-block:: text

   The following databases will be installed:

     pfam Pfam37.4 -> Pfam38.2
       into: /db/pfam/Pfam38.2
       config.yaml key: PFAM_DB

   7 database(s) checked: 5 to install, 1 up to date, 1 skipped.

Options:

- ``databases``: names to update. All of them are checked when none is given.
- ``--yes``: download and install. Without it the command is a dry run.
- ``--no-set-default``: install the releases but leave ``config.yaml`` alone.
  The paths to set are printed at the end.
- ``--timeout``: seconds to wait for each source before giving up.
- The install options of the per-database commands also apply:
  ``--download-runtime``, ``--overwrite``, ``-e/--env_path``, ``-p/--profile``,
  and the resource, Snakemake and SLURM overrides.

Behavior:

- Each release is installed into a new folder beside the configured one, by the
  same workflow ``drakkar database <name>`` runs. Each install keeps its own
  Snakemake working directory, so ``drakkar logging -o <release_dir>``,
  ``drakkar unlock -o <release_dir>``, the failure report and
  ``database_versions.yaml`` work exactly as they do for a single install.
- The installs run sequentially, and one that fails does not stop the others.
  Failed databases are listed at the end with the command that inspects them,
  and ``config.yaml`` is never repointed at a release that failed to install.
- Databases that are up to date are left alone. Databases DRAKKAR does not
  install, such as GTDB, are reported with the newest release available and
  skipped, as are any whose source could not be reached.
- Unless ``--no-set-default`` is given, each installed release becomes the
  default in ``config.yaml``. Existing outputs were built with the earlier
  releases, so rerunning a directory that holds them needs
  ``--allow-database-change`` or the affected outputs deleted. See
  `Cross-run consistency check`_.

Database preflight checks
-------------------------

Every workflow command that uses databases validates them before Snakemake is
launched, so an incomplete or swapped database fails in seconds instead of
after hours of compute.

Installation check
^^^^^^^^^^^^^^^^^^

DRAKKAR verifies that each database the requested module needs is present and
that none of its artifacts are missing or empty. For managed releases this
covers every file the installer produces, including the pressed HMM indices,
the KEGG hierarchy JSON, the KOfam ``ko_list`` cutoff table, and the Pfam EC
mapping table. Only the databases the run actually needs are checked, so
``--annotation-type kegg`` does not require a Pfam release.

If something is missing, DRAKKAR names the exact files and prints the command
that reinstalls the release:

.. code-block:: console

   ERROR: Required databases are missing or incomplete:
     KEGG/KOfam (KEGG_DB): /db/kofams/2026-02-01
       missing or empty: /db/kofams/2026-02-01/kofams_ko_list.tsv
       reinstall with: drakkar database kegg --directory /db/kofams --version 2026-02-01

This matters most for artifacts whose absence would otherwise be silent. A
missing ``ko_list`` stops the KEGG merge with an explicit error, but a missing
KEGG hierarchy JSON would simply drop every EC annotation without any warning.

Use ``--skip-database-check`` to launch without this validation.

Cross-run consistency check
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each run records the databases it used in its ``drakkar_<run_id>.yaml`` run
metadata file:

.. code-block:: yaml

   databases:
     KEGG_DB:
       configured: /db/kofams/2026-02-01
       release: '2026-02-01'
       requested_version: '2026-02-01'
       source_version: kofam archive 2026-02-01

When a later run reuses the same output directory, DRAKKAR compares the
configured databases against the most recent recorded run. Output directories
created before provenance recording fall back to the database paths in
``annotating/annotation_manifest.yaml``; if neither is available, the check is
skipped.

Snakemake profiles use ``rerun-trigger: mtime``, so changing a database in
``config.yaml`` never invalidates existing outputs. Results built with the old
release are silently kept and merged with results built with the new one. To
prevent that, DRAKKAR stops a run whose databases changed **and** whose output
directory still holds files built with the earlier release:

.. code-block:: console

   ERROR: Databases changed since run 20260825-162931, and outputs built with the earlier releases are still present:
     KEGG_DB: 2026-02-01 (/db/kofams/2026-02-01) -> 2026-04-01 (/db/kofams/2026-04-01)
       built with the earlier release: /output/annotating/kegg/PV-171_bin_31.tsv

Resolve it in one of three ways:

- restore the earlier database values in ``config.yaml`` (``drakkar config --edit``)
- delete the listed outputs so they are rebuilt against the new release
- rerun with ``--allow-database-change`` to knowingly mix releases in the
  directory, which is recorded in the run metadata

A database change in a directory that holds no outputs built from it is
reported as information only and does not stop the run.

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

   $ drakkar logging -o drakkar_output --failures

.. code-block:: console

   $ drakkar logging -o drakkar_output --run 20260503-101530 --paths

Options:

- ``-o/--output``: output directory to inspect.
- ``--run``: specific run ID (``YYYYMMDD-HHMMSS``) or
  ``drakkar_<run_id>.yaml`` file name.
- ``--summary``: print only the parsed workflow summary.
- ``--failures``: print the full failure report, including per-job log
  excerpts.
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
- If failures are detected, the failure report described below is printed as
  well.
- If the output directory is locked, run ``drakkar logging -o <output_dir>``
  before using ``drakkar unlock`` or ``--overwrite``.

.. _failure-report:

Failure report
--------------

When Snakemake stops after failures, DRAKKAR prints a tabular failure report
before exiting, and writes the same information to
``drakkar_<run_id>_failures.tsv`` in the root of the output directory, next to
``drakkar_<run_id>.yaml`` and ``drakkar_<run_id>_resources.yaml``. The report
can be printed again at any time with
``drakkar logging -o <output_dir> --failures``.

.. code-block:: text

   ================================ FAILURE REPORT ================================
   Failed jobs: 2 across 2 rule(s) (3 failed attempts in total)
   Failed once but completed after a retry: 1

   Rule                       Target             Att  Reason          Detail
   samtools_stats             PR04534              1  command-error   samtools index: failed to open "PR04534.bam"
   singlem                    PR04533              1  timeout         hit the SLURM time limit (72 min requested)

   =============================== WHAT TO DO NEXT ================================
   [command-error] 1 job(s) in samtools_stats: the tool called by the rule exited with a non-zero status.
     -> Inspect the job log of the failed rule and fix the cause before relaunching.
     -> Example job log: .snakemake/slurm_logs/rule_samtools_stats/PR04534/45360507.log
   [timeout] 1 job(s) in singlem: the job hit its SLURM wall-time limit.
     -> Relaunch the same drakkar command with a larger --time-multiplier (e.g. --time-multiplier 2).

   Verdict: at least one failure needs to be inspected and fixed before relaunching.

Report contents:

- One row per failed job, identified by rule and target (the sample, assembly,
  or MAG taken from the Snakemake wildcards), with the number of failed
  attempts.
- A failure category derived from the SLURM job state, the Snakemake error
  message, and the job log of the failed rule: ``timeout``,
  ``out-of-memory``, ``node-failure``, ``cancelled``, ``storage``,
  ``missing-input``, ``missing-output``, ``incomplete``, ``locked``,
  ``command-error``, or ``unknown``.
- A short detail line: the requested runtime or memory for resource failures,
  and the last error line of the job log for tool failures.
- Jobs that failed but succeeded on a later attempt are counted separately and
  marked as ``recovered`` in the TSV instead of being listed as failures.
- Workflow-level errors that are not tied to a single job, such as
  ``MissingInputException`` or a locked directory, are listed below the table.

The verdict line summarizes whether relaunching is enough:

- ``relaunching the same drakkar command should be enough``: only transient
  failures (node failures, cancellations, missing outputs) were detected.
- ``relaunch with the resource multipliers suggested above``: all failures were
  timeouts or out-of-memory kills, so ``--time-multiplier`` or
  ``--memory-multiplier`` should be raised.
- ``at least one failure needs to be inspected and fixed before relaunching``:
  a tool error, missing input, or storage problem needs attention first.

Completed outputs are preserved, so a relaunch resumes where the workflow
stopped.

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
- ``-a/--annotations``: annotation tables, taxonomy, provenance manifest, and
  QC summary.
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
- ``annotating/``: annotation tables and provenance sidecars; see
  :doc:`annotation_tables` for the gene-table schema and 1.x migration guide.
- ``expressing/``: expression outputs.
- ``dereplicating/``: dereplicated genomes in dereplication-only mode.
- ``benchmark/``: per-SLURM-job resource tables written after each workflow run.
- ``drakkar_<run_id>.yaml``: workflow run metadata.
- ``drakkar_<run_id>_resources.yaml``: root-level SLURM resource-efficiency
  summary for the run (CPU time, memory peaks, and efficiency ratios).
- ``drakkar_<run_id>_failures.tsv``: root-level table of failed jobs, written
  only when a run fails (see :ref:`failure-report`).
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
