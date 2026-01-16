![alt text](drakkar.png "DRAKKAR by the AlberdiLab")

# DRAKKAR

DRAKKAR is a Snakemake-based, modular genome-resolved metagenomics pipeline optimized for the Globe Institute's HPC Mjolnir. It runs on Slurm, supports re-entrant modules, and produces standardized outputs for downstream analysis. Full docs: https://drakkar.readthedocs.io/

## Highlights

- Modular workflows: run the full pipeline or just the parts you need.
- HPC-friendly: Snakemake + Slurm with resource-aware rules.
- Reproducible outputs: consistent folder structure and metadata logs.
- Flexible inputs: directories or sample info tables.

## Quickstart

```
module load drakkar/1.0.0
drakkar complete -f input_info.tsv -o drakkar_output
```

## Workflow overview

DRAKKAR is organized into independent modules. Use `drakkar complete` to chain them.

- **Preprocessing**: quality filtering, optional host removal.
- **Cataloging**: assembly, binning, and coverage mapping.
- **Profiling**: dereplication and abundance estimation.
- **Annotating**: taxonomic and functional annotations.
- **Expressing**: metatranscriptomics mapping to annotated genes.
- **Dereplicating**: dereplication only, no read mapping.
- **Transfer**: SFTP transfer of selected outputs.

## Modules

### Preprocessing

Prepares raw metagenomic reads for downstream modules by filtering, trimming, and (optionally) removing host DNA. It produces cleaned paired reads and summary statistics.

Key options:
- `-i/--input`: input directory containing raw reads.
- `-f/--file`: sample info table (TSV) instead of an input directory.
- `-o/--output`: output directory (defaults to current working directory).
- `-r/--reference`: host reference genome file for host removal.
- `-e/--env_path`: conda environments directory (shared).
- `-p/--profile`: Snakemake profile (default: slurm).

```
drakkar preprocessing {arguments}
```

### Cataloging

Builds assemblies and bins contigs into MAGs. It supports individual assemblies, co-assemblies, and coverage-based mapping strategies for binning.

Key options:
- `-i/--input`: input directory (preprocessed reads or sample info output).
- `-f/--file`: sample info table.
- `-o/--output`: output directory.
- `-m/--mode`: assembly modes (e.g. `individual,all`).
- `-c/--multicoverage`: map coverage groups across assemblies.
- `-e/--env_path`: conda environments directory.
- `-p/--profile`: Snakemake profile.

```
drakkar cataloging {arguments}
```

### Profiling

Quantifies genomes or pangenomes by mapping reads to a dereplicated reference set and summarizing abundance metrics.

Key options:
- `-b/--bins_dir`: directory of MAGs/bins.
- `-B/--bins_file`: text file listing MAG/bin paths.
- `-r/--reads_dir`: directory of metagenomic reads.
- `-R/--reads_file`: sample info table with reads.
- `-o/--output`: output directory.
- `-t/--type`: `genomes` or `pangenomes` (default: genomes).
- `-f/--fraction`: compute microbial fraction with singlem.
- `-a/--ani`: dereplication ANI threshold (default: 0.98).
- `-e/--env_path`: conda environments directory.
- `-p/--profile`: Snakemake profile.

```
drakkar profiling {arguments}
```

### Annotating

Annotates dereplicated MAGs with taxonomic and/or functional labels and produces per-genome tables.

Key options:
- `-b/--bins_dir`: directory of MAGs/bins.
- `-B/--bins_file`: text file listing MAG/bin paths.
- `-o/--output`: output directory.
- `--annotation-type`: `taxonomy`, `function`, or `taxonomy,function`.
- `-e/--env_path`: conda environments directory.
- `-p/--profile`: Snakemake profile.

```
drakkar annotating {arguments}
```

### Expressing

Maps metatranscriptomics reads to annotated genes and reports expression metrics by gene and MAG.

Key options:
- `-b/--bins_dir`: directory of MAGs/bins.
- `-B/--bins_file`: text file listing MAG/bin paths.
- `-r/--reads_dir`: directory of transcriptome reads.
- `-R/--reads_file`: sample info table with transcriptome reads.
- `-o/--output`: output directory.
- `-e/--env_path`: conda environments directory.
- `-p/--profile`: Snakemake profile.

```
drakkar expressing {arguments}
```

### Dereplicating

Runs only the dereplication step from profiling and writes dereplicated MAGs to `dereplicating/final`.

Key options:
- `-b/--bins_dir`: directory of MAGs/bins.
- `-B/--bins_file`: text file listing MAG/bin paths.
- `-o/--output`: output directory.
- `-a/--ani`: dereplication ANI threshold (default: 0.98).
- `-e/--env_path`: conda environments directory.
- `-p/--profile`: Snakemake profile.

```
drakkar dereplicating {arguments}
```

### Transfer

Transfers selected outputs via SFTP while preserving the original folder structure.

Key options:
- `--host` and `--user`: SFTP connection (optional if `--erda` is used).
- `--port`: SFTP port (default: 22).
- `-i/--identity`: SSH private key for authentication.
- `-l/--local-dir`: local output directory to transfer.
- `-r/--remote-dir`: remote destination directory.
- `--erda`: use ERDA defaults (`io.erda.dk` and the default ERDA user).
- `--all`: transfer the entire output directory.
- `--data`: transfer everything except `.snakemake`.
- `--results`: transfer the union of `-a/-m/-p/-b`.
- `-a/--annotations`: annotation outputs.
- `-m/--mags`: dereplicated MAGs.
- `-p/--profile`: profiling outputs.
- `-b/--bins`: cataloging bins (recursive).
- `-v/--verbose`: log each transfer on screen.

```
drakkar transfer {arguments}
```

## Complete mode

Run the full pipeline in one command:

```
drakkar complete {arguments}
```

## Inputs

You can provide inputs as a directory of reads or as a sample info table.

### Directory inputs (minimum)

```
drakkar complete -i {input_path} -o {output_path}
```

### Sample info table

A tab-separated file with optional columns depending on module usage:

|sample|rawreads1|rawreads2|reference_name|reference_path|coassembly|coverage|
|---|---|---|---|---|---|---|
|sample1|path/sample1_1.fq.gz|path/sample1_2.fq.gz|ref1|path/ref1.fna|assembly1,all|coverage1|
|sample2|path/sample2_1.fq.gz|path/sample2_2.fq.gz|ref1|path/ref1.fna|assembly2,all|coverage2|

Notes:
- Input read files can be local paths or remote URLs (http/https/ftp).
- The `coassembly` column defines which samples are pooled for co-assembly.
- `-m individual` adds per-sample assemblies in addition to any co-assemblies.
- If `--multicoverage` is set, samples sharing a coverage group map to each other's assemblies.
- Co-assemblies are not compatible with `--multicoverage`.

## Usage examples

### Run complete with a sample info file

```
drakkar complete -f {info_file} -o {output_path} -m individual
```

### Run preprocessing with a host reference

```
drakkar preprocessing -i {input_path} -o {output_path} -r {genome_path}
```

### Run profiling without microbial fraction

```
drakkar profiling -b {bins_dir} -R {reads_file} -o {output_path}
```

### Run dereplication only

```
drakkar dereplicating -b {bins_dir} -o {output_path}
```

### Transfer results to ERDA

```
drakkar transfer --erda -l {output_path} -r /remote/path --results -v
```

## Outputs

All modules write into the output directory you provide. Key locations:

- `preprocessing/` cleaned reads and summaries.
- `cataloging/` assemblies, bins, and metadata.
- `profiling_genomes/` dereplication, mapping, and abundance tables.
- `annotating/` annotation tables.
- `expressing/` expression outputs.
- `dereplicating/` dereplicated genomes (dereplication-only mode).

Each run also writes a metadata file `drakkar_YYYYMMDD-HHMMSS.yaml` capturing CLI arguments and run context.

## Transfer flags

- `--all`: transfer the entire output directory.
- `--data`: transfer everything except `.snakemake`.
- `--results`: transfer the union of `-a/-m/-p/-b`.
- `-a/--annotations`: annotation tables.
- `-m/--mags`: dereplicated MAGs.
- `-p/--profile`: profiling tables.
- `-b/--bins`: cataloging bins (recursive).
- `--erda`: use `io.erda.dk` with the default ERDA user.
- `-v/--verbose`: log each transfer on screen.

## Support and docs

- Docs: https://drakkar.readthedocs.io/
- Issues: https://github.com/alberdilab/drakkar
