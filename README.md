![DRAKKAR](drakkar.png "DRAKKAR by the AlberdiLab")

# DRAKKAR

[![Documentation Status](https://readthedocs.org/projects/drakkar/badge/?version=latest)](https://drakkar.readthedocs.io/)
[![GitHub release](https://img.shields.io/github/v/release/alberdilab/drakkar)](https://github.com/alberdilab/drakkar/releases)

DRAKKAR is a modular, genome-resolved metagenomics workflow built with Snakemake for HPC environments. It is designed for re-entrant execution, standardized outputs, and end-to-end analysis from raw reads to MAG profiling, annotation, and downstream inspection.

Full documentation lives at https://drakkar.readthedocs.io/. The README is intentionally short; detailed command usage, input formats, database management, outputs, and troubleshooting are documented on Read the Docs.

## What DRAKKAR can do

- Run a complete genome-resolved metagenomics workflow or only the modules you need.
- Preprocess reads with quality filtering, optional host removal, microbial fraction, and Nonpareil.
- Assemble and bin metagenomes into MAG catalogs.
- Dereplicate MAGs and quantify genomes or pangenomes across samples.
- Annotate MAGs taxonomically and functionally with configurable annotation subsets.
- Map metatranscriptomes to annotated genes for expression analysis.
- Install and manage supported annotation databases from the CLI.
- Inspect workflow logs, recover from locked runs, and transfer results by SFTP.

## Quick start

Install from GitHub:

```bash
pip install git+https://github.com/alberdilab/drakkar.git
```

Run the full workflow from a sample table:

```bash
drakkar complete -f input_info.tsv -o drakkar_output
```

Sample tables can point to local read files, remote read URLs, or a paired-end
ENA/SRA run accession in an `accession` column such as `ERR4303216`.

On clusters that provide a module build, load the module first:

```bash
module load drakkar/<version>
```

## Documentation

Start here on Read the Docs:

- Installation: https://drakkar.readthedocs.io/en/latest/installation.html
- User guide and input formats: https://drakkar.readthedocs.io/en/latest/usage.html
- Workflow modules: https://drakkar.readthedocs.io/en/latest/workflows.html
- Operations, status, logging, transfer, and troubleshooting: https://drakkar.readthedocs.io/en/latest/operations.html
- CLI reference: https://drakkar.readthedocs.io/en/latest/api.html

## Main commands

- `drakkar complete`: chain the main workflow modules in one run.
- `drakkar preprocessing`: clean reads and optionally remove host DNA.
- `drakkar cataloging`: assemble, bin, and summarize MAG catalogs.
- `drakkar profiling`: dereplicate genomes and estimate abundance.
- `drakkar annotating`: run taxonomy and selected functional annotations.
- `drakkar expressing`: quantify microbial gene expression.
- `drakkar dereplicating`: run dereplication without read mapping.
- `drakkar database`: install or update supported annotation databases.
- `drakkar status`: show rule and sample progress for a workflow run.
- `drakkar logging`: inspect run metadata and Snakemake logs.
- `drakkar config`: view or edit the installed workflow configuration.
- `drakkar transfer`: transfer selected outputs via SFTP.

## Support

- Documentation: https://drakkar.readthedocs.io/
- Source code and releases: https://github.com/alberdilab/drakkar
- Contact: antton.alberdi@sund.ku.dk
