![alt text](drakkar.png "DRAKKAR by the AlberdiLab")

**Drakkar** is a genome-resolved metagenomics pipeline optimised for running in the Globe Institute's Mjolnir HPC. It is built in a modular fashion, so that the entire workflow or only parts of it can be executed.

## Quickstart

```
module load drakkar/1.0.0
drakkar preprocessing -i spider_data -o drakkar_test -r spider_data/GCA_044660205.1_ASM4466020v1_genomic.fna
drakkar preprocessing -f input_info.tsv -o drakkar_test
```

## Modules

* **Preprocessing**: quality-filters the reads and optionally removes host DNA.
```
drakkar preprocessing
```
* **Cataloging**: assembles and bins the metagenomic reads using multiple strategies.
* **Dereplicating**: dereplicates bins using multiple criteria.
* **Annotating**: taxonomicallt and functionally annotates bins using multiple databases.

## Preprocessing module

- Quality-filtering using fastp
- Reference genome indexing
- Reference genome mapping
- Metagenomic and host genomic data outputting
