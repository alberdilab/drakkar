Installation
============

DRAKKAR can be installed from PyPI or directly from the GitHub repository, or
accessed as an HPC module on supported systems (e.g., Mjolnir). DRAKKAR 2.0
requires Python 3.10 or newer.

The Python package contains the CLI and workflow files; it is not a standalone
installation of the complete bioinformatics stack. Workflow execution targets
module-based HPC systems and requires an Environment Modules/Lmod ``module``
command. Snakemake and other module names, along with database paths, are
defined in ``drakkar/workflow/config.yaml``. Conda environments provide the
dependencies for rules that declare them.

Install From PyPI
-----------------

.. code-block:: console

   $ pip install drakkar

Install From GitHub
-------------------

Use the GitHub URL when you intentionally want the current development version:

.. code-block:: console

   $ pip install git+https://github.com/alberdilab/drakkar.git

If you are working in a managed HPC environment, prefer using a dedicated
Conda environment:

.. code-block:: console

   $ conda create -n drakkar python=3.12 pip
   $ conda activate drakkar
   $ pip install drakkar

Use on Mjolnir (module)
-----------------------

If your cluster provides a DRAKKAR module, load it before running any command:

.. code-block:: console

   $ module load drakkar/<version>

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

Housekeeping: retiring old environments
---------------------------------------

Snakemake names every deployed environment after a hash of its definition, so
each time a ``workflow/envs/*.yaml`` file changes a new environment is created
and the previous one stays behind. Over several DRAKKAR versions a shared
environment directory therefore keeps growing.

List what the directory contains, and which environments the installed DRAKKAR
version no longer uses:

.. code-block:: console

   $ drakkar environments -e /shared/drakkar_envs --list

Each environment is reported as ``in use`` (its definition matches one shipped
by the installed version), ``orphan`` (built from a superseded definition),
``incomplete`` (a failed deployment), or ``unknown`` (a directory that does not
look like a Snakemake environment, which is only reported and never removed).

Delete the reclaimable ones. Without ``--yes`` this is a dry run that only
prints what would be removed:

.. code-block:: console

   $ drakkar environments -e /shared/drakkar_envs --prune
   $ drakkar environments -e /shared/drakkar_envs --prune --yes

Add ``--no-size`` to skip the directory size computation, which can be slow on
shared parallel filesystems.

.. warning::

   Pruning is relative to the DRAKKAR version you run it with. If two DRAKKAR
   versions share one environment directory, pruning from the older one deletes
   the environments the newer one needs.
