Installation
============

DRAKKAR can be installed with pip or accessed as an HPC module on supported
systems (e.g., Mjolnir). The core workflow is executed with Snakemake and uses
Conda environments for tool dependencies.

Install with pip
----------------

.. code-block:: console

   $ pip install drakkar

If you are working in a managed HPC environment, prefer using a dedicated
Conda environment:

.. code-block:: console

   $ conda create -n drakkar python=3.12
   $ conda activate drakkar
   $ pip install drakkar

Use on Mjolnir (module)
-----------------------

If your cluster provides a DRAKKAR module, load it before running any command:

.. code-block:: console

   $ module load drakkar/1.0.0

Verify installation
-------------------

.. code-block:: console

   $ drakkar --help

Optional: pre-create environments
---------------------------------

You can pre-create Conda environments for all workflows (useful on clusters
where environment creation during a run is slow):

.. code-block:: console

   $ drakkar environments --profile slurm

If you want to use a shared environment directory, pass ``-e/--env_path`` to
commands or set it in your workflow config.
