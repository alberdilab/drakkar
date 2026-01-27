CLI Reference
=============

This page summarizes DRAKKAR's command-line interface. For worked examples
and guidance, see :doc:`usage`.

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

Common options (most modules):

- ``-o/--output``: output directory (default: current directory).
- ``-e/--env_path``: shared Conda environment directory.
- ``-p/--profile``: Snakemake profile (default: slurm).

Module-specific options are detailed in :doc:`usage`.
