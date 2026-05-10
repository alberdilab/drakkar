import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from drakkar.cli_context import ERROR, RESET, WORKFLOW_RUN_COMMANDS
from drakkar.output import print
from drakkar.output_paths import validate_launch_metadata_directory

def get_modules_to_run(command):
    if command == "complete":
        return ["preprocessing", "cataloging", "profiling", "annotating"]
    if command:
        return [command]
    return []

def build_snakemake_log_path(output_dir, run_id):
    return Path(output_dir) / "log" / f"drakkar_{run_id}.snakemake.log"

def build_benchmark_paths(output_dir, run_id):
    benchmark_dir = Path(output_dir) / "benchmark"
    base_name = f"drakkar_{run_id}"
    return {
        "dir": benchmark_dir,
        "jobs": benchmark_dir / f"{base_name}.jobs.tsv",
        "rules": benchmark_dir / f"{base_name}.rules.tsv",
        "summary": Path(output_dir) / f"{base_name}_resources.yaml",
    }

def load_metadata_file(metadata_path):
    if not metadata_path:
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except OSError:
        return None

def write_launch_metadata(args, output_dir, env_path=None):
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
    snakemake_log_path = None
    if args.command in WORKFLOW_RUN_COMMANDS:
        snakemake_log_path = build_snakemake_log_path(output_path, run_id)
        benchmark_paths = build_benchmark_paths(output_path, run_id)
        try:
            snakemake_log_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"{ERROR}ERROR:{RESET} Cannot create log directory for Drakkar run metadata: {snakemake_log_path.parent}")
            print(f"{exc.__class__.__name__}: {exc}")
            return None
    else:
        benchmark_paths = None
    metadata = {
        "run_id": run_id,
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
    if snakemake_log_path is not None:
        metadata["snakemake_log"] = str(snakemake_log_path.resolve())
    if benchmark_paths is not None:
        metadata["benchmark_jobs"] = str(benchmark_paths["jobs"].resolve())
        metadata["benchmark_rules"] = str(benchmark_paths["rules"].resolve())
        metadata["benchmark_summary"] = str(benchmark_paths["summary"].resolve())
        metadata["benchmark_status"] = "pending"
    metadata_path = output_path / f"drakkar_{run_id}.yaml"
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
        raise subprocess.CalledProcessError(return_code, command)

    finalize_launch_metadata(run_info, "success", 0, current_workflow=workflow_name)
    if output_dir:
        generate_run_benchmark(output_dir, metadata_path=metadata_path, quiet=True)
    return return_code
