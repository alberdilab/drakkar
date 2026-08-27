"""Preflight validation and cross-run provenance tracking for Drakkar databases.

Two related checks run before a workflow is launched:

* ``check_database_artifacts`` verifies that every database the requested
  module needs is installed and complete. Without it a missing artifact only
  surfaces deep inside the workflow, after hours of compute, and some missing
  artifacts (a KEGG hierarchy JSON, a Pfam EC table) degrade the results
  silently instead of raising.
* ``check_database_provenance`` compares the databases configured for this run
  against the ones recorded by earlier runs in the same output directory.
  Snakemake profiles use ``rerun-trigger: mtime``, so a database swap never
  invalidates existing outputs: results built with two different releases would
  otherwise be merged into one output directory without any trace.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from drakkar.cli_context import ERROR, INFO, RESET, config_vars
from drakkar.database_registry import (
    MANAGED_DATABASES,
    database_artifact_path,
    database_release_from_config,
)
from drakkar.output import print
from drakkar.run_logs import discover_run_metadata

HMM_PRESSED_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")

# Artifacts that must exist inside a managed release directory. These mirror the
# files produced by the prepare_database rules in workflow/rules/databases.smk.
MANAGED_REQUIRED_ARTIFACTS = {
    "kegg": ("", *HMM_PRESSED_SUFFIXES, ".json", "_ko_list.tsv"),
    "cazy": ("", *HMM_PRESSED_SUFFIXES),
    "pfam": ("", *HMM_PRESSED_SUFFIXES, "_ec.tsv"),
    "vfdb": ("", ".dbtype", ".index", ".tsv"),
    "amr": ("", *HMM_PRESSED_SUFFIXES, ".tsv"),
}

GENE_ANNOTATION_OUTPUTS = (
    "annotating/final/*_genes.tsv",
    "annotating/gene_annotations.tsv.xz",
)
CLUSTER_ANNOTATION_OUTPUTS = (
    "annotating/final/*_clusters.tsv",
    "annotating/cluster_annotations.tsv.xz",
)

# Annotation component (as normalized by normalize_annotation_type) -> database.
# "stale_outputs" lists the outputs that were built with the database and must
# be removed before a different release can be used consistently.
ANNOTATION_DATABASE_REQUIREMENTS = {
    "kegg": {
        "config_key": "KEGG_DB",
        "database": "kegg",
        "source": "kegg",
        "label": "KEGG/KOfam",
        "stale_outputs": ("annotating/kegg/*.tsv", *GENE_ANNOTATION_OUTPUTS),
    },
    "cazy": {
        "config_key": "CAZY_DB",
        "database": "cazy",
        "source": "cazy",
        "label": "CAZy/dbCAN HMMs",
        "stale_outputs": ("annotating/cazy/*/dbCAN_hmm_results.tsv", *GENE_ANNOTATION_OUTPUTS),
    },
    "pfam": {
        "config_key": "PFAM_DB",
        "database": "pfam",
        "source": "pfam",
        "label": "Pfam",
        "stale_outputs": ("annotating/pfam/*.tsv", *GENE_ANNOTATION_OUTPUTS),
    },
    "virulence": {
        "config_key": "VFDB_DB",
        "database": "vfdb",
        "source": "vfdb",
        "label": "VFDB",
        "stale_outputs": ("annotating/vfdb/*.txt", *GENE_ANNOTATION_OUTPUTS),
    },
    "amr": {
        "config_key": "AMR_DB",
        "database": "amr",
        "source": "amr",
        "label": "NCBIfam-AMRFinder",
        "stale_outputs": ("annotating/amr/*.tsv", *GENE_ANNOTATION_OUTPUTS),
    },
    "defense": {
        "config_key": "DEFENSEFINDER_DB",
        "database": None,
        "source": "defensefinder",
        "label": "DefenseFinder models",
        "stale_outputs": ("annotating/defensefinder/*", *GENE_ANNOTATION_OUTPUTS),
    },
    "dbcan": {
        "config_key": "DBCAN_DB",
        "database": None,
        "source": "dbcan",
        "label": "dbCAN",
        "stale_outputs": ("annotating/dbcan/*", *CLUSTER_ANNOTATION_OUTPUTS),
    },
    "antismash": {
        "config_key": "ANTISMASH_DB",
        "database": None,
        "source": "antismash",
        "label": "antiSMASH",
        "stale_outputs": ("annotating/antismash/*", *CLUSTER_ANNOTATION_OUTPUTS),
    },
    "mobile": {
        "config_key": "GENOMAD_DB",
        "database": None,
        "source": "genomad",
        "label": "geNomad",
        "stale_outputs": ("annotating/genomad/*", *CLUSTER_ANNOTATION_OUTPUTS),
    },
}

AMR_WORKFLOW_DATABASES = (
    (
        "AMRFINDER_DB",
        "AMRFinderPlus database",
        ("amr/raw/amrfinder/*", "amr/amr_*.tsv.xz", "amr/manifest.yaml"),
        "amrfinderplus",
    ),
    (
        "CARD_DB",
        "CARD database loaded for local RGI use",
        ("amr/raw/rgi/*", "amr/amr_*.tsv.xz", "amr/manifest.yaml"),
        "card",
    ),
    (
        "GENOMAD_DB",
        "geNomad database",
        ("amr/raw/genomad/*", "amr/mobility_regions.tsv.xz", "amr/amr_mobility.tsv.xz", "amr/manifest.yaml"),
        None,
    ),
)

# Annotation manifest source name -> config key, used when an output directory
# predates provenance recording but already holds an annotation manifest.
MANIFEST_SOURCE_TO_CONFIG_KEY = {
    definition["source"]: definition["config_key"]
    for definition in ANNOTATION_DATABASE_REQUIREMENTS.values()
    if definition["source"]
}


@dataclass(frozen=True)
class DatabaseRequirement:
    """One database a module needs, and the outputs already built from it."""

    config_key: str
    label: str
    configured: str
    database: str | None = None
    stale_outputs: tuple[str, ...] = ()


def _configured_value(config, key):
    value = (config or {}).get(key)
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def gtdb_config_key(gtdb_version=None):
    version = str(gtdb_version or "").strip()
    return f"GTDB_DB_{version}" if version else "GTDB_DB"


def _requirement(config, config_key, label, database=None, stale_outputs=()):
    configured = _configured_value(config, config_key)
    if not configured:
        return None
    return DatabaseRequirement(
        config_key=config_key,
        label=label,
        configured=configured,
        database=database,
        stale_outputs=tuple(stale_outputs),
    )


def annotating_requirements(annotation_type, gtdb_version=None, config=None):
    """Databases needed by the annotating module for the requested components."""
    config = config_vars if config is None else config
    components = {item.strip() for item in str(annotation_type or "").split(",") if item.strip()}
    requirements = []
    if "taxonomy" in components:
        requirement = _requirement(
            config,
            gtdb_config_key(gtdb_version),
            "GTDB-Tk reference data",
            stale_outputs=("annotating/gtdbtk/*", "annotating/genome_taxonomy.tsv"),
        )
        if requirement:
            requirements.append(requirement)
    for component, definition in ANNOTATION_DATABASE_REQUIREMENTS.items():
        if component not in components:
            continue
        requirement = _requirement(
            config,
            definition["config_key"],
            definition["label"],
            database=definition["database"],
            stale_outputs=definition["stale_outputs"],
        )
        if requirement:
            requirements.append(requirement)
    return requirements


def amr_requirements(config=None):
    """Databases that are mandatory for the dedicated AMR workflow."""
    config = config_vars if config is None else config
    return [
        DatabaseRequirement(
            config_key=config_key,
            label=label,
            configured=_configured_value(config, config_key) or "",
            database=database,
            stale_outputs=stale_outputs,
        )
        for config_key, label, stale_outputs, database in AMR_WORKFLOW_DATABASES
    ]


def module_requirements(command, args, config=None):
    """Databases needed by the module the user launched."""
    config = config_vars if config is None else config
    requirements = []
    seen = set()

    def add(requirement):
        if requirement and requirement.config_key not in seen:
            seen.add(requirement.config_key)
            requirements.append(requirement)

    fraction = bool(getattr(args, "fraction", False))
    if command in ("preprocessing", "complete") and fraction:
        add(_requirement(config, "SINGLEM_DB", "SingleM metapackage"))
    if command in ("cataloging", "profiling", "dereplicating", "complete"):
        add(_requirement(config, "CHECKM2_DB", "CheckM2 database"))
    if command in ("profiling", "complete") and fraction:
        add(_requirement(config, "SINGLEM_DB", "SingleM metapackage"))
    if command in ("annotating", "complete"):
        for requirement in annotating_requirements(
            getattr(args, "annotation_type", None),
            getattr(args, "gtdb_version", None),
            config=config,
        ):
            add(requirement)
    if command == "amr":
        for requirement in amr_requirements(config=config):
            add(requirement)
    return requirements


def _is_populated(path):
    path = Path(path)
    if path.is_dir():
        return any(path.iterdir())
    return path.is_file() and path.stat().st_size > 0


def missing_artifacts(requirement):
    """Return the artifacts of one database that are absent or empty."""
    if not str(requirement.configured).strip():
        return [f"{requirement.config_key} is not configured in config.yaml"]
    configured = Path(requirement.configured)
    if requirement.database == "amrfinderplus":
        release_dir = database_release_from_config(requirement.database, configured)
        if not release_dir.is_dir():
            return [str(release_dir)]
        required = (
            release_dir / "AMRProt.fa",
            release_dir / "AMRProt.fa.phr",
            release_dir / "AMRProt.fa.pin",
            release_dir / "AMRProt.fa.psq",
            release_dir / "database_format_version.txt",
            release_dir / "version.txt",
            release_dir / "fam.tsv",
        )
        return [str(path) for path in required if not _is_populated(path)]
    if requirement.database in MANAGED_REQUIRED_ARTIFACTS:
        release_dir = database_release_from_config(requirement.database, configured)
        if not release_dir.is_dir():
            return [str(release_dir)]
        return [
            str(artifact)
            for artifact in (
                database_artifact_path(requirement.database, configured, suffix)
                for suffix in MANAGED_REQUIRED_ARTIFACTS[requirement.database]
            )
            if not _is_populated(artifact)
        ]
    if requirement.config_key == "CARD_DB":
        local_database = configured / "localDB"
        if not local_database.is_dir() or not _is_populated(local_database):
            return [str(local_database)]
        return []
    if not configured.exists() or not _is_populated(configured):
        return [str(configured)]
    return []


def reinstall_command(requirement):
    """The drakkar command that reinstalls a managed database release."""
    if requirement.database not in MANAGED_DATABASES:
        return None
    release_dir = database_release_from_config(requirement.database, Path(requirement.configured))
    return (
        f"drakkar database {requirement.database} "
        f"--directory {release_dir.parent} --version {release_dir.name}"
    )


def check_database_artifacts(requirements, skip=False):
    """Report databases that are missing artifacts. Returns True when usable."""
    if skip:
        return True
    problems = [
        (requirement, missing)
        for requirement in requirements
        if (missing := missing_artifacts(requirement))
    ]
    if not problems:
        return True

    print(f"{ERROR}ERROR:{RESET} Required databases are missing or incomplete:")
    for requirement, missing in problems:
        print(f"  {requirement.label} ({requirement.config_key}): {requirement.configured}")
        for path in missing:
            print(f"    missing or empty: {path}")
        command = reinstall_command(requirement)
        if command:
            print(f"    reinstall with: {command}")
    print(f"{INFO}Fix the databases above, or rerun with --skip-database-check to launch anyway.{RESET}")
    return False


def _release_record(configured, database=None):
    """Describe one configured database for the run metadata."""
    if not str(configured).strip():
        return {"configured": "", "release": ""}
    path = Path(configured)
    if database in MANAGED_DATABASES:
        path = database_release_from_config(database, path)
    record = {"configured": str(path), "release": path.name}
    versions_file = path / "database_versions.yaml"
    if not versions_file.is_file():
        versions_file = path.parent / "database_versions.yaml"
    if versions_file.is_file():
        try:
            installed = yaml.safe_load(versions_file.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError):
            installed = {}
        for key in ("requested_version", "source_version"):
            if installed.get(key):
                record[key] = installed[key]
    return record


def collect_database_provenance(requirements):
    """Map each configured database to the record stored in the run metadata."""
    return {
        requirement.config_key: _release_record(requirement.configured, requirement.database)
        for requirement in requirements
    }


def manifest_provenance(output_dir):
    """Read database provenance from an existing annotation manifest."""
    manifest_path = Path(output_dir) / "annotating" / "annotation_manifest.yaml"
    if not manifest_path.is_file():
        return None, {}
    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return None, {}
    records = {}
    for source, described in (manifest.get("databases") or {}).items():
        config_key = MANIFEST_SOURCE_TO_CONFIG_KEY.get(source)
        if not config_key or not isinstance(described, dict):
            continue
        configured = described.get("configured_path")
        if not configured:
            continue
        record = {"configured": str(configured), "release": Path(str(configured)).name}
        for key in ("requested_version", "source_version"):
            if described.get(key):
                record[key] = described[key]
        records[config_key] = record
    if not records:
        return None, {}
    return f"annotation manifest {manifest_path}", records


def previous_database_provenance(output_dir):
    """Find the most recent database provenance recorded in an output directory."""
    for metadata_path, metadata in discover_run_metadata(output_dir):
        databases = metadata.get("databases")
        if databases:
            run_id = metadata.get("run_id") or metadata_path.name
            return f"run {run_id}", databases
    return manifest_provenance(output_dir)


def _record_identity(record):
    if not isinstance(record, dict):
        return None
    return (
        str(record.get("configured") or ""),
        str(record.get("requested_version") or record.get("release") or ""),
    )


def compare_database_provenance(current, previous):
    """Config keys recorded by both runs whose database release differs."""
    changes = []
    for config_key, record in (current or {}).items():
        earlier = (previous or {}).get(config_key)
        if earlier is None:
            continue
        if _record_identity(earlier) != _record_identity(record):
            changes.append((config_key, earlier, record))
    return changes


def stale_outputs_for(output_dir, requirement):
    """Existing outputs in the directory that were built with another release."""
    output_path = Path(output_dir)
    matches = []
    for pattern in requirement.stale_outputs:
        matches.extend(sorted(str(match) for match in output_path.glob(pattern)))
    return matches


def _describe_record(record):
    if not isinstance(record, dict):
        return str(record)
    version = record.get("requested_version") or record.get("release")
    configured = record.get("configured", "")
    return f"{version} ({configured})" if version else str(configured)


def check_database_provenance(output_dir, requirements, current, allow_change=False):
    """Block a run that would mix database releases with existing outputs."""
    source, previous = previous_database_provenance(output_dir)
    if not previous:
        return True
    changes = compare_database_provenance(current, previous)
    if not changes:
        return True

    by_key = {requirement.config_key: requirement for requirement in requirements}
    blocking = []
    for config_key, earlier, record in changes:
        requirement = by_key.get(config_key)
        stale = stale_outputs_for(output_dir, requirement) if requirement else []
        if stale:
            blocking.append((config_key, earlier, record, stale))

    if not blocking:
        print(f"{INFO}INFO:{RESET} Databases changed since {source}, but no outputs were built with the earlier releases:")
        for config_key, earlier, record in changes:
            print(f"  {config_key}: {_describe_record(earlier)} -> {_describe_record(record)}")
        return True

    label = "WARNING" if allow_change else "ERROR"
    colour = INFO if allow_change else ERROR
    print(f"{colour}{label}:{RESET} Databases changed since {source}, and outputs built with the earlier releases are still present:")
    for config_key, earlier, record, stale in blocking:
        print(f"  {config_key}: {_describe_record(earlier)} -> {_describe_record(record)}")
        preview = stale[:5]
        for path in preview:
            print(f"    built with the earlier release: {path}")
        if len(stale) > len(preview):
            print(f"    ... and {len(stale) - len(preview)} more")
    if allow_change:
        print(f"{INFO}Continuing because --allow-database-change was given. The output directory will mix database releases.{RESET}")
        return True

    print("Snakemake reruns on file timestamps only, so these outputs will not be rebuilt for the new release.")
    print("Choose one of:")
    print("  - restore the earlier database values in config.yaml (drakkar config --edit)")
    print("  - delete the outputs listed above so they are rebuilt with the new release")
    print("  - rerun with --allow-database-change to knowingly mix releases in this directory")
    return False
