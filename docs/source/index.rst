DRAKKAR documentation
=====================

**DRAKKAR** is a modular, Snakemake-based workflow for genome-resolved
metagenomics. It is designed for HPC use, re-entrant execution, and consistent
output organization across preprocessing, cataloging, profiling, annotation,
expression, and downstream troubleshooting workflows.

If you are new to DRAKKAR, start with :doc:`installation` and then continue to
the :doc:`usage` guide.

Find what you need
------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - If you want to...
     - Start here
   * - Install DRAKKAR or prepare a working environment
     - See :doc:`installation` for GitHub installation, HPC module use, and
       environment preparation.
   * - Understand the basic input formats and how a run is organized
     - See :doc:`usage` for quickstart commands, sample table structure, core
       concepts, and guide navigation.
   * - Choose and run a workflow module
     - See :doc:`workflows` for the complete workflow and all analysis modules
       such as preprocessing, cataloging, profiling, annotating, expressing,
       dereplicating, and inspecting.
   * - Interpret annotation tables or migrate results from Drakkar 1.x
     - See :doc:`annotation_tables` for examples, schemas, source values,
       migration mappings, and analysis recipes.
   * - Manage databases, logs, config, transfers, outputs, or troubleshooting
     - See :doc:`operations` for operational commands and maintenance tasks.
   * - Summarize a finished run, or query its results as a database
     - See :doc:`report` for the ``drakkar report`` command and the schema of
       the report database.
   * - Check the command-line command list at a glance
     - See :doc:`api` for the CLI reference.

Documentation map
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Start here

   installation
   usage

.. toctree::
   :maxdepth: 2
   :caption: Run workflows

   workflows
   operations

.. toctree::
   :maxdepth: 1
   :caption: Reference

   annotation_tables
   report
   api
