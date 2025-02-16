![alt text](drakkar.png "DRAKKAR by the AlberdiLab")

**Drakkar** is a genome-resolved metagenomics pipeline optimised for running in the Globe Institute's Mjolnir HPC. It is built in a modular fashion, so that the entire workflow or only parts of it can be executed.

## Quickstart

```
module load drakkar/1.0.0
drakkar complete -f input_info.tsv -o drakkar_output
```

## Modules
**DRAKKAR** is a modular software that allows executing each section of the genome-resolved metagenomic pipeline independently.

* **Preprocessing**: quality-filters the reads and optionally removes host DNA.
```
drakkar preprocessing {arguments}
```
* **Cataloging**: assembles and bins the metagenomic reads using multiple strategies.
```
drakkar cataloging {arguments}
```
* **Dereplicating**: dereplicates bins using multiple criteria.
```
drakkar dereplicating {arguments}
```
* **Annotating**: taxonomicallt and functionally annotates bins using multiple databases.
```
drakkar annotating {arguments}
```

## Complete mode
All the modules of **DRAKKAR** can be run together by using the ***complete*** mode.
```
drakkar complete {arguments}
```

## Preprocessing module

- Quality-filtering using fastp
- Reference genome indexing
- Reference genome mapping
- Metagenomic and host genomic data outputting
