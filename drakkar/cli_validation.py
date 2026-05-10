import argparse
import re
from datetime import datetime, timezone

from drakkar.cli_context import (
    CATALOGING_BINNER_ALIASES,
    CATALOGING_BINNER_ORDER,
    DEFAULT_CATALOGING_BINNERS,
    ERROR,
    RESET,
    config_vars,
)
from drakkar.database_registry import MANAGED_DATABASES, normalize_managed_database_name
from drakkar.output import print

def normalize_annotation_type(annotation_type):
    functional_components = {
        "kegg", "cazy", "pfam", "virulence", "amr", "signalp",
        "dbcan", "antismash", "defense", "mobile"
    }
    gene_components = {"kegg", "cazy", "pfam", "virulence", "amr", "signalp"}
    aliases = {"vfdb": "virulence", "genomad": "mobile"}
    allowed = {
        "taxonomy", "function", "genes", "network",
        *functional_components
    }
    option_order = [
        "taxonomy", "function", "genes", "network",
        "kegg", "cazy", "pfam", "virulence", "amr", "signalp",
        "dbcan", "antismash", "defense", "mobile"
    ]
    items = [aliases.get(item.strip().lower(), item.strip().lower()) for item in annotation_type.split(",") if item.strip()]
    invalid = [item for item in items if item not in allowed]
    if not items or invalid:
        print(f"{ERROR}ERROR:{RESET} --annotation-type must be a comma-separated list including taxonomy, function, genes, kegg, cazy, pfam, virulence, amr, signalp, dbcan, antismash, defense, mobile, and/or network.")
        return None

    expanded = set(items)
    if "function" in expanded:
        expanded.update(functional_components)
    if "genes" in expanded:
        expanded.update(gene_components)

    normalized = [opt for opt in option_order if opt in expanded]
    return ",".join(normalized)

def available_gtdb_versions(config=None):
    source = config if config is not None else config_vars
    source = source or {}
    versions = []
    for key, value in source.items():
        match = re.fullmatch(r"GTDB_DB_(\d+)", str(key))
        if match and value:
            versions.append(match.group(1))
    return sorted(set(versions), key=lambda version: int(version), reverse=True)

def validate_gtdb_version(version, config=None):
    if not version:
        return None
    version = str(version).strip()
    if not re.fullmatch(r"\d+", version):
        print(f"{ERROR}ERROR:{RESET} --gtdb-version must be a GTDB release number such as 232.")
        return None

    versions = available_gtdb_versions(config)
    if version not in versions:
        supported = ", ".join(versions) if versions else "none configured"
        print(f"{ERROR}ERROR:{RESET} --gtdb-version {version} is not configured. Available versions: {supported}.")
        return None
    return version

def validate_database_version(version):
    version = (version or "").strip()
    if not version or version in {".", ".."} or "/" in version or "\\" in version:
        print(f"{ERROR}ERROR:{RESET} --version must be a single folder name, not a path.")
        return None
    return version

def validate_download_runtime(value):
    try:
        runtime = int(value)
    except (TypeError, ValueError):
        print(f"{ERROR}ERROR:{RESET} --download-runtime must be a positive integer number of minutes.")
        return None
    if runtime <= 0:
        print(f"{ERROR}ERROR:{RESET} --download-runtime must be a positive integer number of minutes.")
        return None
    return runtime

def normalize_cataloging_binners(binners):
    if not binners:
        return DEFAULT_CATALOGING_BINNERS

    raw_items = [item.strip().lower() for item in str(binners).split(",") if item.strip()]
    if not raw_items:
        print(f"{ERROR}ERROR:{RESET} --binners must include at least one of: {', '.join(CATALOGING_BINNER_ORDER)}.")
        return None
    if "all" in raw_items:
        return DEFAULT_CATALOGING_BINNERS

    invalid = [item for item in raw_items if item not in CATALOGING_BINNER_ALIASES]
    if invalid:
        print(
            f"{ERROR}ERROR:{RESET} --binners contains unsupported value(s): {', '.join(invalid)}. "
            f"Options are: {', '.join(CATALOGING_BINNER_ORDER)}."
        )
        return None

    selected = {CATALOGING_BINNER_ALIASES[item] for item in raw_items}
    return ",".join(binner for binner in CATALOGING_BINNER_ORDER if binner in selected)

def positive_int(value):
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed

def add_resource_multiplier_arguments(parser):
    parser.add_argument(
        "--memory-multiplier",
        type=positive_int,
        default=1,
        help=(
            "Multiply Snakemake memory requests before applying SNAKEMAKE_MAX_GB. "
            "Default: 1"
        ),
    )
    parser.add_argument(
        "--time-multiplier",
        type=positive_int,
        default=1,
        help=(
            "Multiply Snakemake runtime requests before applying SNAKEMAKE_MAX_TIME. "
            "Default: 1"
        ),
    )

def add_benchmark_argument(parser):
    parser.add_argument(
        "--skip-benchmark",
        action="store_true",
        help="Skip generating post-run SLURM resource benchmark outputs",
    )

def resource_config(memory_multiplier=1, time_multiplier=1):
    return f"memory_multiplier={memory_multiplier} time_multiplier={time_multiplier} "

def default_resource_args(memory_multiplier=1, time_multiplier=1):
    max_mem_mb = int(config_vars.get("SNAKEMAKE_MAX_GB", 1024)) * 1024
    max_runtime = int(config_vars.get("SNAKEMAKE_MAX_TIME", 14 * 24 * 60))
    default_mem_mb = min(max_mem_mb, 8 * 1024 * memory_multiplier)
    default_runtime = min(max_runtime, 10 * time_multiplier)
    return f"--default-resources mem_mb={default_mem_mb} runtime={default_runtime} "

def default_database_version(database_name):
    if database_name == "vfdb":
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")
    return None

def validate_managed_database_version(database_name, version):
    if not version:
        default_version = default_database_version(database_name)
        if default_version:
            return default_version
        print(f"{ERROR}ERROR:{RESET} --version is required for {database_name}.")
        return None
    version = validate_database_version(version)
    if not version:
        return None
    if database_name == "kegg":
        try:
            parsed = datetime.strptime(version, "%Y-%m-%d")
        except ValueError:
            print(f"{ERROR}ERROR:{RESET} KEGG --version must be an archive date in YYYY-MM-DD format, e.g. 2026-02-01")
            return None
        return parsed.strftime("%Y-%m-%d")
    return version
