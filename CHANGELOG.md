# Changelog

This project tracks release notes here from this point forward.

## [Unreleased]

### Added

- No unreleased changes yet.

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
