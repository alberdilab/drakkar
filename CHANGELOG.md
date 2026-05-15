# Changelog

This project tracks release notes here from this point forward.

## [Unreleased]

### Added

- No unreleased changes yet.

## [1.6.9] - 2026-05-15

### Changed

- GTDB-Tk taxonomy annotation now resolves the configured GTDB database at rule execution time, reports the selected config key and path, and uses the configured GTDB-Tk module for classification.

### Fixed

- ENA/SRA accession FASTQ downloads now validate cached and downloaded files against ENA byte counts when available, retry truncated files, and reject clearly imbalanced paired FASTQ sizes before Snakemake starts when ENA byte counts are unavailable.

## [1.6.8] - 2026-05-10

### Added

- Added `drakkar status` to show latest or selected run progress by rule and sample, with optional run-directory/YAML selection, rule-only/sample-only views, and `--complete` helper-rule detail.

## [1.6.7] - 2026-05-10

### Changed

- Increased cataloging runtime scaling for SemiBin2 and ComeBin so their requested time is based on full assembly size instead of half assembly size.

### Fixed

- Fixed SLURM benchmark capture for current Snakemake executor logs such as `Job <id> has been submitted with SLURM jobid <id>`, and normalized embedded `sbatch` responses before querying `sacct`.

## [1.6.6] - 2026-05-10

### Added

- Added stricter input-transfer preflight checks before Snakemake starts: URL downloads now retry up to five times, empty cached/downloaded files are rejected, SFTP URLs are supported via `curl`, and missing or empty input files/directories are reported together with a clear stop message.

### Changed

- Split the large CLI and utility modules into focused implementation modules while keeping the existing `drakkar.cli` and `drakkar.utils` import surfaces compatible.

## [1.6.5] - 2026-05-09

### Added

- Added `comebin` as a cataloging binner and added `-b/--binners` to `drakkar cataloging` and `drakkar complete` to select comma-separated binners from `metabat`, `maxbin`, `semibin`, and `comebin`; the default remains all binners.

## [1.6.4] - 2026-05-07

### Fixed

- When no assembly column is present in the sample info file (or all assembly values are empty) and no `--mode` flag is given, DRAKKAR now defaults to individual assemblies instead of failing with an empty `assembly_to_samples.json`.

## [1.6.3] - 2026-05-07

### Changed

- Cataloging now prioritises preprocessed reads when resolving input paths from a sample info file: explicit `preprocessedreads1`/`preprocessedreads2` columns are used first, then reads found in `preprocessing/final/` from a prior preprocessing run, and raw `rawreads1`/`rawreads2` or `accession` paths are only used as a last resort. This allows the same input file to carry coverage/assembly grouping metadata while still assembling and mapping against quality-filtered reads.

## [1.6.2] - 2026-05-05

### Added

- Added pruned GTDB-Tk output trees for taxonomy annotation so ``annotating/bacteria.tree`` and optional ``annotating/archaea.tree`` keep only the input MAG genomes and exclude GTDB reference tips.

### Changed

- ENA/SRA and URL file downloads now retry up to three times with exponential backoff (5 s, 15 s, 45 s) before giving up on a single file.
- All file-downloading loops (samples, references, bins, MAGs, preprocessed reads, transcriptomes) now continue past individual download failures and report the full list of failed files at the end instead of stopping on the first error.

## [1.6.1] - 2026-05-04

### Changed

- Made SLURM resource benchmarking explicitly default-on across workflow commands, added a shared ``--skip-benchmark`` opt-out flag, and ensured empty benchmark tables are still written when no submitted jobs are detected so the benchmark output is visible by default.
- Broadened Snakemake SLURM submission parsing to recognize additional submission log styles such as ``Submitted batch job ...`` when building run benchmarks.
- Allowed sample tables for preprocessing, complete, profiling, and expressing to use an ``accession`` column with paired-end ENA/SRA run accessions, which DRAKKAR now resolves and downloads into the existing forward/reverse read caches automatically.
- Styled the screen-session startup hint so the inline ``screen -S mysession`` command is visually marked as code in Rich output and still remains clear in plain-text output.

## [1.6.0] - 2026-05-04

### Added

- Added a first SLURM-backed benchmark layer that parses submitted non-local Snakemake jobs, infers relaunch attempts, joins them to `sacct`, writes per-run benchmark tables under `benchmark/`, emits a root-level `drakkar_<run_id>_resources.yaml` summary for every workflow run, and exposes a compact resource-efficiency summary through `drakkar logging`.
- Added a richer root-level `profiling_genomes.tsv` summary with input reads, input bases, mapped reads, mapping percentage, and mapped bases per sample, and added a new root-level `dereplicating.tsv` summary for profiling and dereplicating runs with input/output bin counts plus average completeness and contamination when quality metadata are available.

## [1.5.4] - 2026-05-04

### Changed

- Animated the DRAKKAR banner sequence in interactive terminals so the ship appears first, the logo follows after a short pause, and the remaining intro or help content appears after a second pause.
- Allowed ``-B/--bins_file`` in profiling, dereplicating, annotating, inspecting, and expressing to point to a remote manifest URL, downloading the manifest and any referenced remote MAG/bin FASTA files into the local cache.
- Made ``--quality`` genome-name matching more flexible so MAG/bin names are accepted with or without FASTA extensions and still map back to the canonical genome filenames used internally.
- Added ``drakkar update --skip-deps`` to let users skip reinstalling Python dependencies during an update when they only want to refresh the package itself.

## [1.5.3] - 2026-05-03

### Changed

- Styled the update-success ASCII banner with the same ship color as the main DRAKKAR banner and rendered the ship version badge with a separate accent color when Rich output is available.

### Fixed

- Made the cataloging Binette checkpoint tolerate empty ``temporary_files/diamond_result.tsv`` or ``diamond_result.tsv.gz`` files by exporting a header-only ``final_bins_quality_reports.tsv`` instead of failing the workflow.

## [1.5.2] - 2026-05-03

### Changed

- Stopped ``drakkar logging`` from dumping the Snakemake log by default, and replaced that section with a short guide showing how to inspect summaries, excerpts, paths, full logs, and available runs.
- Added an explicit ``--excerpt`` mode to ``drakkar logging`` so failure excerpts or fallback tails are shown only on request.

### Fixed

- Ensured cataloging still writes ``cataloging/final/all_bin_metadata.csv`` when Binette yields zero bins, by materializing an empty metadata table with the expected columns.

## [1.5.1] - 2026-05-03

### Added

- Added a plain-help fallback note explaining when Rich-styled help is unavailable because the ``rich`` dependency is missing.

### Changed

- Rendered the DRAKKAR ship and logo above CLI help output, consolidated the top-level workflow table into ``Data Generation and Analysis``, and expanded workflow descriptions in the top-level help menu.

### Fixed

- Changed ``drakkar update`` to reinstall package dependencies instead of skipping them, so upgraded environments can pick up required runtime packages such as ``rich``.
- Normalized release-changelog generation so published versions always remain separated by a blank line.

## [1.5.0] - 2026-05-03

### Added

- Added parsed execution summaries to ``drakkar logging``, including planned jobs, observed rule executions, workflow progress, and detected error types, plus a ``--summary`` mode for summary-only output.
- Added a version-aware success banner at the end of ``drakkar update`` showing the version that was just installed.

### Changed

- Moved DRAKKAR ASCII art assets into ``drakkar/ascii.py`` so display helpers in ``utils.py`` only handle formatting and printing.
- Reworked ``drakkar --help`` and subcommand help into grouped Rich layouts with workflow categories, command families, option sections, and concrete examples.

## [1.4.0] - 2026-05-03

### Added

- Accepted ``assembly`` as the preferred sample-table column for cataloging assembly groups, while keeping legacy ``coassembly`` support.
- Added ``drakkar logging`` plus persistent per-run Snakemake logs under ``log/drakkar_<run_id>.snakemake.log`` to help troubleshoot failed or locked workflow directories.

## [1.3.5] - 2026-05-02

### Added

- Added `--memory-multiplier` and `--time-multiplier` to Snakemake-running Drakkar commands.
- Added `SNAKEMAKE_MAX_TIME` to cap runtime requests, defaulting to 20160 minutes (14 days).

### Changed

- Routed explicit workflow memory and runtime resources through capped scaling helpers, and scaled default Snakemake resources for rules without explicit resources.

## [1.3.4] - 2026-05-02

### Fixed

- Changed the SemiBin2 cataloging table to read reclustered `output_bins/*.fa` files instead of pre-reclustering `contig_bins.tsv`.
- Counted MetaBAT2, MaxBin2, and SemiBin2 bins in `cataloging.tsv` from Binette input quality reports instead of contig-to-bin table rows.

## [1.3.3] - 2026-05-01

### Changed

- Rendered argparse help with Rich panels and tables for command and subcommand help.

### Fixed

- Forced Drakkar CLI truecolor output for normal terminal streams so inherited non-color environment settings do not hide the Rich theme.
- Merged SingleM microbial fraction rows whose sample names include read-mate suffixes such as `_1` back into the main preprocessing summary rows.

## [1.3.2] - 2026-05-01

### Changed

- Moved the startup version badge to the ship banner, centered the intro against the DRAKKAR logo, and enriched CLI status/help styling through the shared Rich output layer.

## [1.3.1] - 2026-05-01

### Changed

- Styled the startup banner with a logo-integrated version badge, centered intro frame, and Rich colors inspired by Viking sea, iron, and gold tones.

## [1.3.0] - 2026-05-01

### Added

- Cataloging now runs QUAST assembly statistics and samtools flagstat mapping summaries, then writes a root-level `cataloging.tsv` summary table.

### Changed

- Host/metagenomic read and base count sidecar files are now temporary preprocessing intermediates; the retained summary is `preprocessing.tsv`.

## [1.2.3] - 2026-05-01

### Added

- All SingleM `pipe` and microbial fraction steps now use the configured `SINGLEM_DB` metapackage.

## [1.2.2] - 2026-04-30

### Added

- Added `--fraction` to `drakkar preprocessing` so standalone preprocessing can run SingleM and SingleM microbial fraction outputs.
- Added `--nonpareil` to `drakkar preprocessing` and `drakkar complete` to estimate Nonpareil coverage/diversity from preprocessed metagenomic reads.
- Expanded `preprocessing.tsv` with fastp, host/metagenomic, SingleM fraction, and Nonpareil summary columns.

## [1.2.1] - 2026-04-30

### Fixed

- Allowed preprocessing `-r/--reference` and `-x/--reference-index` to accept `http`, `https`, and `ftp` URLs.

## [1.2.0] - 2026-04-30

### Added

- Added `-x/--reference-index` for preprocessing and complete workflows to use tarballs containing a host FASTA plus prebuilt Bowtie2 index files, including tarballs listed in sample table `reference_path` values.

### Fixed

- Added a preflight check for current/output directory access before writing Drakkar run metadata, so protected directories produce a clear CLI error instead of a Python `PermissionError`.

## [1.1.3] - 2026-04-27

### Added

- Added a `.gitignore` for Python caches, test artifacts, build outputs, documentation builds, Snakemake state, and local OS/editor files.
- Added `--gtdb-version` to `drakkar annotating` and `drakkar complete` to select numbered `GTDB_DB_<version>` config entries for taxonomy annotation.

### Fixed

- Capped dynamic Snakemake `mem_mb` resource requests using configurable `SNAKEMAKE_MAX_GB`, defaulting to 1024 GB.

## [1.1.2] - 2026-04-25

### Added

- Added Rich-backed CLI output, prompts, help rendering, and `drakkar --version` for the 1.1.2 release.

## [1.1.1] - 2026-04-24

### Added

- Added `drakkar config --view` and `drakkar config --edit` to inspect or modify the installed `workflow/config.yaml` from the CLI.

### Changed

- Added `--overwrite` to output-writing workflows so a locked output directory can be deleted and recreated after a broken Snakemake session, with an interactive prompt when the flag is not provided.
- Changed `drakkar database amr` to download versioned NCBI AMRFinder releases such as `2025-07-16.1` instead of the floating `latest` endpoint.
- Changed `drakkar database kegg` to download KOfam profiles from versioned KEGG monthly archive directories using `YYYY-MM-DD` archive dates, then merge and `hmmpress` the extracted HMMs.
- Changed the default database download/preparation runtime to `120` minutes for the KEGG archive rule, and added `--download-runtime` to `drakkar database` so this limit can be overridden from the CLI.
- Changed managed annotation database config entries to point to the installed release directory, with the workflow resolving the expected internal files automatically.
- Changed `drakkar database pfam` to download `Pfam-A.hmm.gz` from versioned Pfam release directories such as `Pfam37.4`, then unzip and `hmmpress` the requested release.
- Changed `drakkar database vfdb` so `--version` can be omitted and then defaults to the UTC download date, which is used as the release folder and logged version.
- Fixed `drakkar database cazy` to download the requested upstream dbCAN release from the versioned `Databases/<version>/` endpoint instead of silently saving an HTML landing page.

## [1.1.0] - 2026-04-24

### Added

- Added a new `drakkar database <name>` maintenance command for installing one managed annotation database release at a time (`kegg`, `cazy`, `pfam`, `vfdb`, `amr`).
- Added per-run database version logging in `database_versions.yaml`, including source URLs, requested versions, release directories, timestamps, and SHA256 checksums of installed assets.

### Changed

- Managed database installs now use explicit `--directory` and `--version` arguments, and `--set-default` updates the concrete path stored in `config.yaml` only after a successful installation.
- Restored direct-path config usage for non-managed and externally sourced databases such as `GTDB_DB`.

### Documentation

- Documented the new `database` workflow and its version logging behavior in the README and usage guide.

## [1.0.1] - 2026-04-21

### Changed

- Added release automation scaffolding, including this changelog, package-version wiring, and release helper tooling.
