import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from drakkar import __version__
from drakkar.cli_context import ERROR, RESET, WORKFLOW_RUN_COMMANDS
from drakkar.output import print
from drakkar.output_paths import validate_launch_metadata_directory

def get_modules_to_run(command):
    if command == "complete":
        return ["preprocessing", "cataloging", "profiling", "annotating"]
    if command:
        return [command]
    return []

# Everything a workflow run leaves behind about itself — its metadata, its
# Snakemake log, its failure table and its benchmark roll-up — shares one
# directory, named after the `drakkar logging` command that reads them and
# after the gerund convention the workflow output directories already follow.
# The benchmark artefacts nest one level deeper: only SLURM runs produce them,
# and keeping the roll-up with the tables it summarizes leaves `drakkar_*.yaml`
# in the logging directory matching run metadata and nothing else.
LOGGING_DIRNAME = "logging"
BENCHMARK_DIRNAME = "benchmark"

# Runs launched before the logging directory existed wrote their metadata,
# failure table and benchmark roll-up flat into the output root, their
# Snakemake log into `log/` and their benchmark tables into `benchmark/`, with
# `_resources`/`_failures` joined by an underscore rather than a dot. Those
# directories are still read, so every path builder takes a `legacy` flag and
# every discovery helper searches both layouts.
LEGACY_LOG_DIRNAME = "log"


def logging_dir(output_dir):
    """Return the directory holding one output directory's run artefacts."""
    return Path(output_dir) / LOGGING_DIRNAME

def benchmark_dir(output_dir, legacy=False):
    """Return the directory holding the per-run benchmark artefacts."""
    if legacy:
        return Path(output_dir) / BENCHMARK_DIRNAME
    return logging_dir(output_dir) / BENCHMARK_DIRNAME

def build_metadata_path(output_dir, run_id, legacy=False):
    name = f"drakkar_{run_id}.yaml"
    if legacy:
        return Path(output_dir) / name
    return logging_dir(output_dir) / name

def build_snakemake_log_path(output_dir, run_id, legacy=False):
    name = f"drakkar_{run_id}.snakemake.log"
    if legacy:
        return Path(output_dir) / LEGACY_LOG_DIRNAME / name
    return logging_dir(output_dir) / name

def build_benchmark_paths(output_dir, run_id, legacy=False):
    directory = benchmark_dir(output_dir, legacy=legacy)
    base_name = f"drakkar_{run_id}"
    summary_name = f"{base_name}_resources.yaml" if legacy else f"{base_name}.resources.yaml"
    return {
        "dir": directory,
        "jobs": directory / f"{base_name}.jobs.tsv",
        "rules": directory / f"{base_name}.rules.tsv",
        "summary": (Path(output_dir) / summary_name) if legacy else (directory / summary_name),
    }

def uses_legacy_layout(output_dir, metadata_path):
    """True when a run's metadata sits outside the logging directory.

    A run's artefacts follow its metadata file, so regenerating a benchmark or
    a failure report for a run that predates the logging directory rewrites the
    files beside the ones it already has instead of scattering the set across
    both layouts.
    """
    if not metadata_path:
        return False
    try:
        return Path(metadata_path).resolve().parent != logging_dir(output_dir).resolve()
    except OSError:
        return True

def find_snakemake_log(output_dir, run_id):
    """Return a run's Snakemake log from either layout, or None if absent."""
    for legacy in (False, True):
        candidate = build_snakemake_log_path(output_dir, run_id, legacy=legacy)
        if candidate.exists():
            return candidate
    return None

def load_metadata_file(metadata_path):
    if not metadata_path:
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except OSError:
        return None

def write_launch_metadata(args, output_dir, env_path=None, databases=None):
    output_path = Path(output_dir)
    if not validate_launch_metadata_directory(output_path):
        return None
    try:
        output_path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"{ERROR}ERROR:{RESET} Cannot create output directory for Drakkar run metadata: {output_path}")
        print(f"{exc.__class__.__name__}: {exc}")
        return None
    timestamp = datetime.now(timezone.utc)
    run_id = timestamp.strftime("%Y%m%d-%H%M%S")
    run_logging_dir = logging_dir(output_path)
    try:
        run_logging_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"{ERROR}ERROR:{RESET} Cannot create logging directory for Drakkar run metadata: {run_logging_dir}")
        print(f"{exc.__class__.__name__}: {exc}")
        return None
    snakemake_log_path = None
    if args.command in WORKFLOW_RUN_COMMANDS:
        snakemake_log_path = build_snakemake_log_path(output_path, run_id)
        benchmark_paths = build_benchmark_paths(output_path, run_id)
    else:
        benchmark_paths = None
    metadata = {
        "run_id": run_id,
        "drakkar_version": __version__,
        "timestamp": timestamp.isoformat(),
        "started_at": timestamp.isoformat(),
        "command": args.command,
        "modules": get_modules_to_run(args.command),
        "working_directory": str(Path.cwd()),
        "output_directory": str(output_path.resolve()),
        "arguments": vars(args),
        "argv": sys.argv,
        "status": "prepared",
    }
    if env_path is not None:
        metadata["env_path"] = env_path
    if databases:
        metadata["databases"] = databases
    if snakemake_log_path is not None:
        metadata["snakemake_log"] = str(snakemake_log_path.resolve())
    if benchmark_paths is not None:
        metadata["benchmark_jobs"] = str(benchmark_paths["jobs"].resolve())
        metadata["benchmark_rules"] = str(benchmark_paths["rules"].resolve())
        metadata["benchmark_summary"] = str(benchmark_paths["summary"].resolve())
        metadata["benchmark_status"] = "pending"
    metadata_path = build_metadata_path(output_path, run_id)
    try:
        with open(metadata_path, "w") as f:
            yaml.safe_dump(metadata, f, sort_keys=False)
    except OSError as exc:
        print(f"{ERROR}ERROR:{RESET} Cannot write Drakkar run metadata: {metadata_path}")
        print("Run drakkar from a writable directory or pass -o/--output to a writable output directory.")
        print(f"{exc.__class__.__name__}: {exc}")
        return None
    return {
        "run_id": run_id,
        "metadata_path": metadata_path,
        "snakemake_log_path": snakemake_log_path,
        "benchmark_paths": benchmark_paths,
    }

def update_launch_metadata(metadata_path, **updates):
    if not metadata_path:
        return None
    metadata_path = Path(metadata_path)
    if not metadata_path.exists():
        return None
    metadata = load_metadata_file(metadata_path)
    if metadata is None:
        return None
    metadata.update(updates)
    try:
        with open(metadata_path, "w", encoding="utf-8") as handle:
            yaml.safe_dump(metadata, handle, sort_keys=False)
    except OSError:
        return None
    return metadata

def finalize_launch_metadata(run_info, status, exit_code=None, current_workflow=None):
    if not run_info:
        return None
    payload = {
        "status": status,
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    if exit_code is not None:
        payload["exit_code"] = exit_code
    if current_workflow is not None:
        payload["current_workflow"] = current_workflow
    return update_launch_metadata(run_info["metadata_path"], **payload)

def run_subprocess_with_logging(command, run_info=None, workflow_name=None):
    from drakkar.benchmark import generate_run_benchmark
    from drakkar.failures import report_run_failures

    metadata_path = run_info["metadata_path"] if run_info else None
    log_path = Path(run_info["snakemake_log_path"]) if run_info and run_info.get("snakemake_log_path") else None
    output_dir = None
    if run_info and metadata_path:
        metadata = load_metadata_file(metadata_path) or {}
        output_dir = metadata.get("output_directory")
    if metadata_path:
        update_launch_metadata(
            metadata_path,
            status="running",
            current_workflow=workflow_name,
        )

    log_handle = None
    try:
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_handle = open(log_path, "a", encoding="utf-8")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        if process.stdout is not None:
            try:
                for line in process.stdout:
                    print(line, end="")
                    if log_handle is not None:
                        log_handle.write(line)
            finally:
                process.stdout.close()
        return_code = process.wait()
    except Exception:
        if log_handle is not None:
            log_handle.flush()
            log_handle.close()
        finalize_launch_metadata(run_info, "failed", current_workflow=workflow_name)
        raise
    finally:
        if log_handle is not None and not log_handle.closed:
            log_handle.flush()
            log_handle.close()

    if return_code != 0:
        finalize_launch_metadata(run_info, "failed", return_code, current_workflow=workflow_name)
        if output_dir:
            generate_run_benchmark(output_dir, metadata_path=metadata_path, quiet=True)
            report_run_failures(output_dir, metadata_path=metadata_path)
        raise subprocess.CalledProcessError(return_code, command)

    finalize_launch_metadata(run_info, "success", 0, current_workflow=workflow_name)
    if output_dir:
        generate_run_benchmark(output_dir, metadata_path=metadata_path, quiet=True)
    return return_code
