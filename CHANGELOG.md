# Changelog

This project tracks release notes here from this point forward.

## [Unreleased]

### Added

- No unreleased changes yet.

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
