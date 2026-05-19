CLI Reference
=============

This page summarizes DRAKKAR's command-line interface. For worked examples and
workflow guidance, see :doc:`usage`, :doc:`workflows`, and :doc:`operations`.

Global
------

.. code-block:: console

   $ drakkar <command> [options]

Commands
--------

complete
^^^^^^^^
Run the full pipeline.

preprocessing
^^^^^^^^^^^^^
Quality filtering and optional host removal.

cataloging
^^^^^^^^^^
Assembly, binning, and bin metadata.

profiling
^^^^^^^^^
Dereplication and abundance profiling.

annotating
^^^^^^^^^^
Taxonomic and functional annotation.

expressing
^^^^^^^^^^
Metatranscriptome mapping to annotated genes.

dereplicating
^^^^^^^^^^^^^
Dereplication only (no read mapping).

inspecting
^^^^^^^^^^
Microdiversity inspection workflows.

environments
^^^^^^^^^^^^
Pre-create Conda environments for modules.

unlock
^^^^^^
Unlock a Snakemake working directory.

update
^^^^^^
Reinstall DRAKKAR from the Git repository.

transfer
^^^^^^^^
Transfer outputs via SFTP while preserving folder structure.

Selected options
----------------

Common options (all workflow modules):

- ``-o/--output``: output directory (default: current directory).
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile (default: ``slurm``).
- ``--overwrite``: delete a locked output directory and rerun from scratch.
- ``--skip-benchmark``: skip SLURM resource benchmark collection after the run.
- ``--memory-multiplier N``: scale per-rule memory requests before the cap.
- ``--time-multiplier N``: scale per-rule runtime requests before the cap.

Snakemake override options (all workflow modules):

- ``--snakemake-jobs N``: max concurrent SLURM jobs.
- ``--snakemake-cores N``: max local CPU cores (local executor only).
- ``--snakemake-executor EXECUTOR``: Snakemake executor (e.g. ``slurm``, ``local``).
- ``--snakemake-latency-wait N``: seconds to wait for output files.
- ``--snakemake-retries N``: retry failed jobs N times.
- ``--snakemake-rerun-incomplete``: rerun jobs with incomplete output files.
- ``--snakemake-keep-going``: continue after failures.

SLURM override options (all workflow modules):

- ``--slurm-partition NAME``: SLURM partition/queue.
- ``--slurm-account NAME``: SLURM billing account.
- ``--slurm-constraint EXPR``: node constraint expression.
- ``--slurm-nodes N``: nodes per SLURM job.
- ``--slurm-nodelist NODES``: restrict to specific node(s).
- ``--slurm-extra ARGS``: extra ``sbatch`` arguments passed verbatim.

Module-specific options are detailed in :doc:`workflows`. Full Snakemake and
SLURM management documentation is in :doc:`operations`.
