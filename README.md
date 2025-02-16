# DRAKKAR

██████████   ███████████     █████████   █████   ████ █████   ████   █████████   ███████████

░███░░░░███ ░░███░░░░░███   ███░░░░░███ ░░███   ███░ ░░███   ███░   ███░░░░░███ ░░███░░░░░███

░███   ░░███ ░███    ░███  ░███    ░███  ░███  ███    ░███  ███    ░███    ░███  ░███    ░███

░███    ░███ ░██████████   ░███████████  ░███████     ░███████     ░███████████  ░██████████

░███    ░███ ░███░░░░░███  ░███░░░░░███  ░███░░███    ░███░░███    ░███░░░░░███  ░███░░░░░███

░███    ███  ░███    ░███  ░███    ░███  ░███ ░░███   ░███ ░░███   ░███    ░███  ░███    ░███

██████████   █████   █████ █████   █████ █████ ░░████ █████ ░░████ █████   █████ █████   █████

░░░░░░░░░░   ░░░░░   ░░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░ ░░░░░   ░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░░


**Drakkar** is a genome-resolved metagenomics pipeline optimised for running in the Globe Institute's Mjolnir HPC. It is built in a modular fashion, so that the entire workflow or only parts of it can be executed.

```
conda activate drakkar
conda install anaconda::pip
pip uninstall drakkar -y
pip install git+https://github.com/alberdilab/drakkar.git
drakkar --help
```

```
drakkar preprocessing -i spider_data -o drakkar_test -r spider_data/GCA_044660205.1_ASM4466020v1_genomic.fna
drakkar assembly -i drakkar_test/preprocessing/final -o drakkar_test -m individual
```

## Preprocessing module

- Quality-filtering using fastp
- Reference genome indexing
- Reference genome mapping
- Metagenomic and host genomic data outputting
