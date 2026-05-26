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

def nonnegative_float(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("must be a non-negative number") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be a non-negative number")
    return parsed

def percent_float(value):
    parsed = nonnegative_float(value)
    if parsed > 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
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

def add_snakemake_override_arguments(parser):
    snakemake_group = parser.add_argument_group("Snakemake overrides")
    snakemake_group.add_argument(
        "--snakemake-latency-wait",
        type=positive_int,
        dest="snakemake_latency_wait",
        metavar="N",
        help="Seconds to wait for missing output files before failing. Overrides the profile value (slurm default: 300, local default: 60).",
    )
    snakemake_group.add_argument(
        "--snakemake-jobs",
        type=positive_int,
        dest="snakemake_jobs",
        metavar="N",
        help="Maximum number of concurrent SLURM jobs. Overrides the profile value (default: 100).",
    )
    snakemake_group.add_argument(
        "--snakemake-cores",
        type=positive_int,
        dest="snakemake_cores",
        metavar="N",
        help="Maximum local CPU cores (local executor only). Overrides the profile value.",
    )
    snakemake_group.add_argument(
        "--snakemake-executor",
        dest="snakemake_executor",
        metavar="EXECUTOR",
        help="Snakemake executor to use (e.g. slurm, local). Overrides the profile value.",
    )
    snakemake_group.add_argument(
        "--snakemake-retries",
        type=positive_int,
        dest="snakemake_retries",
        metavar="N",
        help="Number of times to retry a failed job. Overrides the profile value (slurm default: 3).",
    )
    snakemake_group.add_argument(
        "--snakemake-rerun-incomplete",
        dest="snakemake_rerun_incomplete",
        action="store_true",
        default=None,
        help="Force rerun of incomplete jobs from a previous interrupted run.",
    )
    snakemake_group.add_argument(
        "--snakemake-keep-going",
        dest="snakemake_keep_going",
        action="store_true",
        default=None,
        help="Continue executing independent jobs after a failure instead of stopping immediately.",
    )
    slurm_group = parser.add_argument_group("SLURM overrides")
    slurm_group.add_argument(
        "--slurm-partition",
        dest="slurm_partition",
        metavar="NAME",
        help="SLURM partition (queue) to submit jobs to.",
    )
    slurm_group.add_argument(
        "--slurm-account",
        dest="slurm_account",
        metavar="NAME",
        help="SLURM account for job billing.",
    )
    slurm_group.add_argument(
        "--slurm-constraint",
        dest="slurm_constraint",
        metavar="EXPR",
        help="SLURM node constraint expression (e.g. gpu or skylake).",
    )
    slurm_group.add_argument(
        "--slurm-nodes",
        type=positive_int,
        dest="slurm_nodes",
        metavar="N",
        help="Number of nodes per SLURM job. Default: 1.",
    )
    slurm_group.add_argument(
        "--slurm-nodelist",
        dest="slurm_nodelist",
        metavar="NODES",
        help="Specific node or node list to run jobs on (e.g. node01 or node[01-03]).",
    )
    slurm_group.add_argument(
        "--slurm-qos",
        dest="slurm_qos",
        metavar="NAME",
        help="SLURM Quality of Service (QOS) string passed to sbatch --qos.",
    )
    slurm_group.add_argument(
        "--slurm-extra",
        dest="slurm_extra",
        metavar="ARGS",
        help="Arbitrary extra sbatch arguments passed verbatim (e.g. '--mail-type=END --mail-user=you@example.com').",
    )

def build_snakemake_flags(args):
    """Build extra top-level Snakemake CLI flags from override arguments."""
    parts = []
    if getattr(args, "snakemake_latency_wait", None) is not None:
        parts.append(f"--latency-wait {args.snakemake_latency_wait}")
    if getattr(args, "snakemake_jobs", None) is not None:
        parts.append(f"--jobs {args.snakemake_jobs}")
    if getattr(args, "snakemake_cores", None) is not None:
        parts.append(f"--cores {args.snakemake_cores}")
    if getattr(args, "snakemake_executor", None):
        parts.append(f"--executor {args.snakemake_executor}")
    if getattr(args, "snakemake_retries", None) is not None:
        parts.append(f"--retries {args.snakemake_retries}")
    if getattr(args, "snakemake_rerun_incomplete", None):
        parts.append("--rerun-incomplete")
    if getattr(args, "snakemake_keep_going", None):
        parts.append("--keep-going")
    return (" ".join(parts) + " ") if parts else ""

def build_slurm_resource_overrides(args):
    """Build SLURM resource override string to append into --default-resources."""
    parts = []
    if getattr(args, "slurm_partition", None):
        parts.append(f"slurm_partition={args.slurm_partition}")
    if getattr(args, "slurm_account", None):
        parts.append(f"slurm_account={args.slurm_account}")
    if getattr(args, "slurm_constraint", None):
        parts.append(f"slurm_constraint={args.slurm_constraint}")
    if getattr(args, "slurm_nodes", None) is not None:
        parts.append(f"slurm_nodes={args.slurm_nodes}")
    if getattr(args, "slurm_qos", None):
        parts.append(f"slurm_qos={args.slurm_qos}")
    nodelist = getattr(args, "slurm_nodelist", None)
    extra = getattr(args, "slurm_extra", None)
    combined_extra_parts = []
    if nodelist:
        combined_extra_parts.append(f"--nodelist={nodelist}")
    if extra:
        combined_extra_parts.append(extra)
    if combined_extra_parts:
        combined = " ".join(combined_extra_parts)
        parts.append(f"slurm_extra='{combined}'")
    return " ".join(parts)

def add_benchmark_argument(parser):
    parser.add_argument(
        "--skip-benchmark",
        action="store_true",
        help="Skip generating post-run SLURM resource benchmark outputs",
    )

def resource_config(memory_multiplier=1, time_multiplier=1):
    return f"memory_multiplier={memory_multiplier} time_multiplier={time_multiplier} "

def default_resource_args(memory_multiplier=1, time_multiplier=1, slurm_resources=""):
    max_mem_mb = int(config_vars.get("SNAKEMAKE_MAX_GB", 1024)) * 1024
    max_runtime = int(config_vars.get("SNAKEMAKE_MAX_TIME", 14 * 24 * 60))
    default_mem_mb = min(max_mem_mb, 8 * 1024 * memory_multiplier)
    default_runtime = min(max_runtime, 10 * time_multiplier)
    extra = f" {slurm_resources}" if slurm_resources else ""
    return f"--default-resources mem_mb={default_mem_mb} runtime={default_runtime}{extra} "

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
