![alt text](drakkar.png "DRAKKAR by the AlberdiLab")

**DRAKKAR** is a snakemake-based genome-resolved metagenomics pipeline optimised for the Globe Institute's HPC Mjolnir. Snakemake works along with Slurm to run long pipelines using optimal memory and time resources. It is built in a modular fashion, so that the entire workflow or only parts of it can be executed. Extended usage tutorial can be found in https://drakkar.readthedocs.io/

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
* **Annotating**: annotates the bins taxonomically and/or functionally, and/or generates community-scale metabolic networks.
```
drakkar annotating {arguments}
```
* **Profiling**: conducts genome- or pangenome-based quantitative analyses.
```
drakkar profiling {arguments}
```

* **Expressing**: conducts gene expression analysis.
```
drakkar profiling {arguments}
```

## Complete mode
All the modules of **DRAKKAR** can be run together by using the ***complete*** mode.
```
drakkar complete {arguments}
```

## Usage examples

### Without sample info file

#### Minimum usage

```
drakkar complete -i {input_path} -o {output_path}
```

**-i**: path to the folder where the metagenomic sequencing reads are stored.
**-o**: path in which the DRAKKAR outputs will be stored.
Metagenomic reads are not mapped to a host genome, individual assemblies are performed, and genome-based profiling is conducted.

#### With reference genome

```
drakkar complete -i {input_path} -o {output_path} -r {genome_path}
```

**-i**: path to the folder where the metagenomic sequencing reads are stored.
**-o**: path in which the DRAKKAR outputs will be stored.
**-r**: path to the reference genome.
Metagenomic reads are mapped to the host genome individual assemblies are performed, and genome-based profiling is conducted.

#### With reference genome and assembly mode

```
drakkar complete -i {input_path} -o {output_path} -r {genome_path} -m individual,all -t genomes,pangenomes
```

**-i**: path to the folder where the metagenomic sequencing reads are stored.
**-o**: path in which the DRAKKAR outputs will be stored.
**-r**: path to the reference genome.
**-m**: comma-separated list of assembly modes
Metagenomic reads are mapped to the host genome, individual assemblies as well as a single coassembly including all samples are performed, and both genome- and pangenome-based profiling is conducted.

### With sample info file

|sample|rawreads1|rawreads2|reference_name|reference_path|assembly|
|---|---|---|---|---|---|
|sample1|path/sample1_1.fq.gz|path/sample1_2.fq.gz|ref1|path/ref1.fna|assembly1,all|
|sample1|path/sample1_1.fq.gz|path/sample1_2.fq.gz|ref1|path/ref1.fna|assembly1,all|
|sample2|path/sample2_1.fq.gz|path/sample2_2.fq.gz|ref1|path/ref1.fna|assembly2,all|
|sample3|path/sample3_1.fq.gz|path/sample3_2.fq.gz|ref2|path/ref2.fna|assembly2,all|
|sample4|path/sample4_1.fq.gz|path/sample4_2.fq.gz|ref2|path/ref2.fna|assembly2,all|
|sample4|path/sample4_1.fq.gz|path/sample4_2.fq.gz|ref2|path/ref2.fna|assembly2,all|

#### Minimum usage

```
drakkar complete -f {info_file} -o {output_path}
```
All the required information is extracted from the sample info file.

#### Minimum usage

```
drakkar complete -f {info_file} -o {output_path} -m individual
```
Individual assemblies are also conducted on top of the assemblies specified in the sample info file.


## DRAKKAR modules

### Preprocessing

- Quality-filtering using fastp
- Reference genome indexing
- Reference genome mapping
- Metagenomic and host genomic data outputting

### Cataloging

Documentation to be added.

### Profiling

Documentation to be added.
