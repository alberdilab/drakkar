import argparse
import os
import sys
import subprocess
import yaml
import re
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import PurePosixPath
from drakkar.utils import *

###
# Define and read config file
###

PACKAGE_DIR = Path(__file__).parent
CONFIG_PATH = PACKAGE_DIR / "workflow" / "config.yaml"

def load_config():
    """Load fixed variables from config.yaml."""
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH, "r") as f:
            return yaml.safe_load(f)
    return {}

config_vars = load_config()

###
# Define text colors
###

HEADER1 = "\033[1;95m"
ERROR = "\033[1;31m"
INFO = "\033[1;34m"
RESET = "\033[0m"

###
# Define helper functions
###

def normalize_annotation_type(annotation_type):
    allowed = {"taxonomy", "function"}
    items = [item.strip() for item in annotation_type.split(",") if item.strip()]
    invalid = [item for item in items if item not in allowed]
    if not items or invalid:
        print(f"{ERROR}ERROR:{RESET} --annotation-type must be 'taxonomy', 'function', or a comma-separated list of those options (e.g. taxonomy,function).")
        return None
    return ",".join(items)

def validate_path(path_value, label, expect_dir=False):
    if not path_value:
        return True
    if expect_dir:
        exists = os.path.isdir(path_value)
    else:
        exists = os.path.exists(path_value)
    if not exists:
        kind = "directory" if expect_dir else "path"
        print(f"{ERROR}ERROR:{RESET} {label} {kind} not found: {path_value}")
        return False
    return True

def get_modules_to_run(command):
    if command == "complete":
        return ["preprocessing", "cataloging", "profiling", "annotating"]
    if command:
        return [command]
    return []

def write_launch_metadata(args, output_dir, env_path=None):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc)
    metadata = {
        "timestamp": timestamp.isoformat(),
        "command": args.command,
        "modules": get_modules_to_run(args.command),
        "working_directory": str(Path.cwd()),
        "output_directory": str(output_path.resolve()),
        "arguments": vars(args),
        "argv": sys.argv,
    }
    if env_path is not None:
        metadata["env_path"] = env_path
    filename_timestamp = timestamp.strftime("%Y%m%d-%H%M%S")
    metadata_path = output_path / f"drakkar_{filename_timestamp}.yaml"
    with open(metadata_path, "w") as f:
        yaml.safe_dump(metadata, f, sort_keys=False)

def list_files_recursive(base_dir, exclude_snakemake=False):
    for root, dirs, files in os.walk(base_dir):
        if exclude_snakemake:
            dirs[:] = [d for d in dirs if d != ".snakemake"]
        for filename in files:
            yield Path(root) / filename

def collect_transfer_files(base_dir, args):
    selected_files = set()
    missing_paths = []
    warnings = []

    def add_annotations():
        add_file(Path("annotating/cluster_annotations.tsv.xz"))
        add_file(Path("annotating/gene_annotations.tsv.xz"))
        add_file(Path("annotating/genome_taxonomy.tsv"))

    def add_mags():
        add_dir(Path("profiling_genomes/drep/dereplicated_genomes"))

    def add_profile():
        add_file(Path("profiling_genomes/final/bases.tsv"))
        add_file(Path("profiling_genomes/final/counts.tsv"))
        add_file(Path("profiling_genomes.tsv"))

    def add_bins():
        add_dir(Path("cataloging/final"))

    def add_file(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_file():
            selected_files.add(local_path)
        else:
            missing_paths.append(str(local_path))

    def add_dir(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_dir():
            for file_path in local_path.rglob("*"):
                if file_path.is_file():
                    selected_files.add(file_path)
        else:
            missing_paths.append(str(local_path))

    if args.all:
        for file_path in list_files_recursive(base_dir):
            if file_path.is_file():
                selected_files.add(file_path)
        return selected_files, missing_paths, warnings

    if args.data:
        for file_path in list_files_recursive(base_dir, exclude_snakemake=True):
            if file_path.is_file():
                selected_files.add(file_path)

    if args.results:
        add_annotations()
        add_mags()
        add_profile()
        add_bins()

    if args.annotations:
        add_annotations()

    if args.mags:
        add_mags()

    if args.profile:
        add_profile()

    if args.bins:
        add_bins()

    return selected_files, missing_paths, warnings

def build_sftp_batch_commands(files, base_dir, remote_dir):
    remote_root = PurePosixPath(remote_dir)
    mkdirs = set()
    put_cmds = []

    for file_path in sorted(files):
        rel_path = file_path.relative_to(base_dir)
        remote_path = remote_root / PurePosixPath(rel_path.as_posix())
        parent = remote_path.parent
        while parent != remote_root and parent not in mkdirs:
            mkdirs.add(parent)
            parent = parent.parent
        put_cmds.append((file_path, remote_path))

    commands = []
    for directory in sorted(mkdirs, key=lambda p: len(p.parts)):
        if directory == remote_root:
            continue
        commands.append(f'-mkdir "{directory.as_posix()}"')
    for local_path, remote_path in put_cmds:
        commands.append(f'put "{local_path}" "{remote_path.as_posix()}"')

    return "\n".join(commands) + "\n", put_cmds

def run_sftp_transfer(args):
    base_dir = Path(args.local_dir).resolve()
    if args.erda:
        args.host = "io.erda.dk"
        args.user = "antton.alberdi@snm.ku.dk"
    if not args.host or not args.user:
        print(f"{ERROR}ERROR:{RESET} --host and --user are required unless --erda is set.")
        return
    selected_files, missing_paths, warnings = collect_transfer_files(base_dir, args)

    if warnings:
        for warning in warnings:
            print(f"{INFO}INFO:{RESET} {warning}")

    if missing_paths:
        for missing in missing_paths:
            print(f"{ERROR}MISSING:{RESET} {missing}")

    if not selected_files:
        print(f"{ERROR}ERROR:{RESET} No files selected for transfer.")
        return

    batch_commands, put_cmds = build_sftp_batch_commands(selected_files, base_dir, args.remote_dir)
    if args.verbose:
        for local_path, remote_path in put_cmds:
            print(f"{INFO}PUT:{RESET} {local_path} -> {remote_path.as_posix()}")
    sftp_cmd = ["sftp"]
    if args.port:
        sftp_cmd += ["-P", str(args.port)]
    if args.identity:
        sftp_cmd += ["-i", args.identity]
    sftp_cmd += ["-b", "-", f"{args.user}@{args.host}"]

    result = subprocess.run(sftp_cmd, input=batch_commands, text=True, capture_output=True)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        mkdir_only_errors = False
        if stderr:
            error_lines = [line for line in stderr.splitlines() if line.strip()]
            mkdir_only_errors = all("Couldn't create directory" in line for line in error_lines)
        if mkdir_only_errors:
            print(f"{INFO}INFO:{RESET} sftp reported existing directories; continuing.")
        else:
            if stderr:
                print(stderr, file=sys.stderr)
            print(f"{ERROR}ERROR:{RESET} sftp transfer failed with exit code {result.returncode}")

###
# Define workflow launching functions
###

def run_unlock(workflow, output_dir, profile):

    unlock_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked")
    print(f"You can now rerun a new workflow using any drakkar command:")

def run_snakemake_environments(workflow, profile):
    cmd = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {Path.cwd()} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} workflow={workflow} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(cmd, shell=False, check=True)

def run_snakemake_preprocessing(workflow, project_name, output_dir, reference, env_path, profile):

    """ Run the preprocessing workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} reference={reference} "
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
        f"--slurm-delete-logfiles-older-than 0"
    ]

    subprocess.run(snakemake_command, shell=False, check=True)


def run_snakemake_cataloging(workflow, project_name, output_dir, env_path, profile):

    """ Run the cataloging workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

#Screen output control
def run_snakemake_cataloging2(workflow, project_name, output_dir, env_path, profile):

    """ Run the cataloging workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
    ]

    process = subprocess.Popen(snakemake_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    for line in process.stdout:
        print(line, end="")  # Print each line as it comes

    # Wait for the process to finish
    process.wait()

    if process.returncode != 0:
        # Read stderr and print it
        error_message = process.stderr.read()
        if "LockException" in error_message:
            display_unlock()
        else:
            print(f"\nERROR: Snakemake failed with exit code {process.returncode}!", file=sys.stderr)
            sys.exit(1)
    else:
        display_end()

def run_snakemake_profiling(workflow, project_name, profiling_type, output_dir, env_path, profile, fraction, ani):
    """ Run the profiling workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} profiling_type={profiling_type} output_dir={output_dir} fraction={fraction} DREP_ANI={ani} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_dereplicating(workflow, project_name, output_dir, env_path, profile, ani):
    """ Run the dereplicating workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} DREP_ANI={ani} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_annotating(workflow, project_name, annotating_type, output_dir, env_path, profile):
    """ Run the profiling workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} annotating_type={annotating_type} output_dir={output_dir} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_inspecting(workflow, project_name, output_dir, env_path, profile):
    """ Run the profiling workflow """

    snakemake_command = [
        "/bin/bash", "-c", 
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_expressing(workflow, project_name, output_dir, env_path, profile):
    """ Run the expressing workflow """

    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

###
# Main function to launch workflows
###

def main():
    parser = argparse.ArgumentParser(
        description="Drakkar: A Snakemake-based workflow for sequencing analysis",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", help="Available workflows")

    # Define subcommands for each workflow
    subparser_complete = subparsers.add_parser("complete", help="Run the complete workflow")
    subparser_complete.add_argument("-i", "--input", required=False, help="Input directory")
    subparser_complete.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_complete.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_complete.add_argument("-r", "--reference", required=False, help="Reference host genome")
    subparser_complete.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_complete.add_argument("-t", "--type", required=False, default="genomes", help="Either genomes or pangenomes profiling type. Default: genomes")
    subparser_complete.add_argument("--annotation-type", dest="annotation_type", required=False, default="taxonomy,function", help="Taxonomic and/or functional annotations (comma-separated). Default: taxonomy,function")
    subparser_complete.add_argument("-c", "--multicoverage", action="store_true", help="Map samples sharing the same coverage group to each other's individual assemblies")
    subparser_complete.add_argument("--fraction", required=False, action='store_true', help="Calculate microbial fraction using singlem")
    subparser_complete.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_complete.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_complete.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_preprocessing = subparsers.add_parser("preprocessing", help="Run the preprocessing workflow (quality-filtering and host removal)")
    subparser_preprocessing.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_preprocessing.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_preprocessing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_preprocessing.add_argument("-r", "--reference", required=False, help="Reference host genome file (fna.gz)")
    subparser_preprocessing.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_preprocessing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_cataloging = subparsers.add_parser("cataloging", help="Run the cataloging (assembly and binning) workflow")
    subparser_cataloging.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_cataloging.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_cataloging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_cataloging.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_cataloging.add_argument("-c", "--multicoverage", action="store_true", help="Map samples sharing the same coverage group to each other's individual assemblies")
    subparser_cataloging.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_cataloging.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_profiling = subparsers.add_parser("profiling", help="Run the profiling workflow")
    subparser_profiling.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_profiling.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_profiling.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_profiling.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_profiling.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_profiling.add_argument("-t", "--type", required=False, default="genomes", help="Either genomes or pangenomes profiling type. Default: genomes")
    subparser_profiling.add_argument("-f", "--fraction", required=False, action='store_true', help="Calculate microbial fraction using singlem")
    subparser_profiling.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_profiling.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_profiling.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_dereplicating = subparsers.add_parser("dereplicating", help="Run dereplication only (no mapping)")
    subparser_dereplicating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_dereplicating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_dereplicating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_dereplicating.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_dereplicating.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_dereplicating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_annotating = subparsers.add_parser("annotating", help="Run the annotating workflow")
    subparser_annotating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa, .fna or .fasta, optionally including .gz) are stored")
    subparser_annotating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_annotating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_annotating.add_argument("--annotation-type", dest="annotation_type", required=False, default="taxonomy,function", help="Taxonomic and/or functional annotations (comma-separated). Default: taxonomy,function")
    subparser_annotating.add_argument("-e", "--env_path", type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_annotating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_inspecting = subparsers.add_parser("inspecting", help="Run the inspecting workflow")
    subparser_inspecting.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_inspecting.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_inspecting.add_argument("-m", "--mapping_dir", required=False, help="Directory containing the mapping (.bam) files")
    subparser_inspecting.add_argument("-c", "--cov_file", required=False, help="Tab-separated file containing mapping coverage per genome per smaple")
    subparser_inspecting.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_inspecting.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_inspecting.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_expressing = subparsers.add_parser("expressing", help="Run the microbial gene expression workflow")
    subparser_expressing.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_expressing.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_expressing.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_expressing.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_expressing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_expressing.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_expressing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_environments = subparsers.add_parser("environments", help="Pre-create conda environments")
    subparser_environments.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_environments.add_argument("--profile", default="local", choices=["local", "slurm"])

    subparser_unlock = subparsers.add_parser("unlock", help="Unlock snakemake")
    subparser_unlock.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_unlock.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_unlock.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_update = subparsers.add_parser("update", help="Reinstall Drakkar from the Git repo (forces reinstall in this environment)")
    subparser_update.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")

    subparser_transfer = subparsers.add_parser("transfer", help="Transfer outputs via sftp")
    subparser_transfer.add_argument("--host", help="SFTP host")
    subparser_transfer.add_argument("--user", help="SFTP user")
    subparser_transfer.add_argument("--port", type=int, default=22, help="SFTP port")
    subparser_transfer.add_argument("-i", "--identity", help="Path to SSH private key")
    subparser_transfer.add_argument("-l", "--local-dir", required=False, default=os.getcwd(), help="Local output directory. Default is the directory from which drakkar is called.")
    subparser_transfer.add_argument("-r", "--remote-dir", required=True, help="Remote base directory for transfer")
    subparser_transfer.add_argument("--erda", action="store_true", help="Use ERDA SFTP defaults")
    subparser_transfer.add_argument("--all", action="store_true", help="Transfer the entire output folder")
    subparser_transfer.add_argument("--data", action="store_true", help="Transfer the output folder excluding .snakemake")
    subparser_transfer.add_argument("--results", action="store_true", help="Transfer configured results files")
    subparser_transfer.add_argument("-a", "--annotations", action="store_true", help="Transfer annotation outputs")
    subparser_transfer.add_argument("-m", "--mags", action="store_true", help="Transfer dereplicated MAGs")
    subparser_transfer.add_argument("-p", "--profile", action="store_true", help="Transfer profiling outputs")
    subparser_transfer.add_argument("-b", "--bins", action="store_true", help="Transfer cataloging bins")
    subparser_transfer.add_argument("-v", "--verbose", action="store_true", help="Log each transfer on screen")

    args = parser.parse_args()

    # Display ASCII logo before running any command or showing help
    display_drakkar()

    # Check screen session
    check_screen_session()

    path_checks = [
        (getattr(args, "input", None), "Input", True),
        (getattr(args, "file", None), "Sample detail file", False),
        (getattr(args, "bins_dir", None), "Bins directory", True),
        (getattr(args, "bins_file", None), "Bins file", False),
        (getattr(args, "reads_dir", None), "Reads directory", True),
        (getattr(args, "reads_file", None), "Reads file", False),
        (getattr(args, "reference", None), "Reference", False),
        (getattr(args, "cov_file", None), "Coverage file", False),
    ]
    for path_value, label, expect_dir in path_checks:
        if not validate_path(path_value, label, expect_dir):
            return

    ###
    # Declare environment directory
    ###

    env_path = getattr(args, "env_path", None)
    if env_path:
        env_path = env_path
    else:
        if args.command != "update":
            env_path = config_vars['ENVIRONMENTS_DIR']

    if args.command in ("annotating", "complete"):
        normalized_annotation_type = normalize_annotation_type(args.annotation_type)
        if not normalized_annotation_type:
            return
        args.annotation_type = normalized_annotation_type

    if args.command == "transfer":
        output_dir = getattr(args, "local_dir", os.getcwd())
    else:
        output_dir = getattr(args, "output", os.getcwd())
    write_launch_metadata(args, output_dir, env_path=locals().get("env_path"))

    ###
    # Unlock, update or create environments
    ###

    if args.command == "unlock":
        print(f"{HEADER1}UNLOCKING DRAKKAR DIRECTORY...{RESET}", flush=True)
        print(f"", flush=True)
        run_unlock(args.command, args.output, args.profile)

    elif args.command == "transfer":
        print(f"{HEADER1}TRANSFERRING DRAKKAR OUTPUTS...{RESET}", flush=True)
        print(f"", flush=True)
        run_sftp_transfer(args)
        return

    elif args.command == "update":
        pip_cmd = [
                    sys.executable, "-m", "pip", "install",
                    "--upgrade", "--force-reinstall", "--no-deps",
                    "git+https://github.com/alberdilab/drakkar.git"
                ]
        print(f"{HEADER1}UPDATING DRAKKAR...{RESET}", flush=True)
        print(f"", flush=True)
        try:
            update_code = subprocess.run(pip_cmd)
        except Exception as e:
            print(f"Update failed: {e}", file=sys.stderr, flush=True)
            sys.exit(1)

    elif args.command == "environments":
        print(f"{HEADER1}CREATING CONDA ENVIRONMENTS...{RESET}", flush=True)
        print(f"", flush=True)
        run_snakemake_environments(args.command, args.env_path, args.profile)

    else:            
        project_name = os.path.basename(os.path.normpath(args.output))

        # Check if directory is locked
        if is_snakemake_locked(args.output):
            print(f"{ERROR}Error: directory is locked by Snakemake{RESET}", file=sys.stderr)
            print(f"This usually happens when a workflow is interrupted suddenly,")
            print(f"often while Slurm jobs are still running. Make sure all pending")
            print(f"jobs are finished or killed before unlocking the working directory.")
            print(f"{INFO}Run 'drakkar unlock' to unlock the working directory{RESET}")
            sys.exit(2)

    ###
    # Preprocessing
    ###

    if args.command in ("preprocessing", "complete"):
        print(f"", flush=True)
        print(f"{HEADER1}STARTING PREPROCESSING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        # Generate raw data dictionaries

        if args.file and args.input:
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.file and not args.input:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.input and not args.file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_samples_to_json(args.input,args.output)

        else:
            print(f"Please provide either an input directory (-i) or a sample info file (-f)")
            return

        # Generate reference genome dictionaries

        if args.file and args.reference:
            if check_reference_columns(args.file):
                print(f"")
                print(f"Both sample info file and reference genome file were provided.")
                print(f"DRAKKAR will continue with the information provided in the sample info file.")
                file_references_to_json(args.file,args.output)
                REFERENCE = True
            else:
                argument_references_to_json(args.reference,f"{args.output}/data/sample_to_reads1.json",args.output)
                REFERENCE = True

        elif args.file and not args.reference:
            if check_reference_columns(args.file):
                print(f"")
                print(f"DRAKKAR will extract the reference genome information from the sample info file.")
                file_references_to_json(args.file,args.output)
                REFERENCE = True
            else:
                print(f"")
                print(f"No reference genome information was provided in the info file.")
                print(f"DRAKKAR will run without mapping against a reference genome.")
                REFERENCE = False

        elif args.reference and not args.file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reference genome file.")
            argument_references_to_json(args.reference,f"{args.output}/data/sample_to_reads1.json",args.output)
            REFERENCE = True

        else:
            print(f"")
            print(f"Running DRAKKAR without mapping against a reference genome")
            REFERENCE = False

        run_snakemake_preprocessing("preprocessing", project_name, Path(args.output).resolve(), REFERENCE, env_path, args.profile)

    ###
    # Cataloging
    ###

    if args.command in ("cataloging", "complete"):
        print(f"", flush=True)
        print(f"{HEADER1}STARTING CATALOGING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        # Generate cataloging data dictionaries

        if args.file and args.input:
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        elif args.file and not args.input:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_preprocessed_to_json(args.file,args.output)
        elif args.input and not args.file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the preprocessing data")
            if any(os.scandir(f"{args.output}/preprocessing/final")):
                argument_preprocessed_to_json(f"{args.output}/preprocessing/final",args.output)
            else:
                print(f"ERROR: No preprocessed data was found in the output directory.")
                print(f"    Please, provide an input directory or sample info file to proceed.")
                return

        # Generate the assembly dictionary

        with open(f"{args.output}/data/preprocessed_to_reads1.json", "r") as f:
            PREPROCESSED_TO_READS1 = json.load(f)

        samples = list(PREPROCESSED_TO_READS1.keys())
        INDIVIDUAL_MODE=False
        ALL_MODE=False

        if args.mode:
            if "individual" in args.mode:
                INDIVIDUAL_MODE=True
            if "all" in args.mode:
                ALL_MODE=True
        else:
            if not args.file:
                print(f"")
                print(f"No assembly mode (-m) or sample info file (-f) has been provided.")
                print(f"    In consequence, DRAKKAR will run individual assemblies.")
                INDIVIDUAL_MODE=True

        file_assemblies_to_json(args.file,samples,INDIVIDUAL_MODE,ALL_MODE,args.output)

        if args.multicoverage:
            with open(f"{args.output}/data/assembly_to_samples.json", "r") as f:
                assembly_to_samples = json.load(f)
            if any(len(group_samples) > 1 for group_samples in assembly_to_samples.values()):
                print(f"{ERROR}ERROR:{RESET} --multicoverage is not compatible with co-assemblies.")
                print("Co-assemblies already use coverage from the samples used to build the assembly.")
                return
            has_coverage = file_coverages_to_json(args.file, samples, args.output)
            if has_coverage:
                print("Coverage groups loaded from the sample info file.")
            else:
                print("No coverage information provided; all samples will be mapped against all assemblies.")

        run_snakemake_cataloging("cataloging", project_name, Path(args.output).resolve(), env_path, args.profile)

    ###
    # Profiling
    ###

    if args.command in ("profiling", "complete"):
        print(f"", flush=True)
        print(f"{HEADER1}STARTING PROFILING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        # Prepare bin dictionaries

        if args.command in ("profiling"):
            bins_dir = args.bins_dir
            bins_file = args.bins_file

        if args.command in ("complete"):
            bins_dir = None
            bins_file = None

        if bins_dir and bins_file:
            print(f"Both bin path file and bin directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the bin info file.")
            print(bins_file)
            file_bins_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"DRAKKAR will run with the bins in the bin directory.")
            path_bins_to_json(bins_dir,args.output)
        else:
            print(f"No bin information was provided. DRAKKAR will try to guess the location of the bins.")
            if os.path.exists(f"{args.output}/cataloging/final/all_bin_paths.txt"):
                file_bins_to_json(f"{args.output}/cataloging/final/all_bin_paths.txt",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an bin file (-B) or directory (-b).")
                return

        # Prepare read dictionaries

        if args.command in ("profiling"):
            reads_dir = args.reads_dir
            reads_file = args.reads_file

        if args.command in ("complete"):
            reads_dir = None
            reads_file = None

        if reads_dir and reads_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_preprocessed_to_json(reads_file,args.output)
        elif reads_file and not reads_dir:
            print(f"")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_preprocessed_to_json(reads_file,args.output)
        elif reads_dir and not reads_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reads from the input directory.")
            argument_preprocessed_to_json(reads_dir,args.output)
        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the reads.")
            if any(os.scandir(f"{args.output}/preprocessing/final")):
                argument_preprocessed_to_json(f"{args.output}/preprocessing/final",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        run_snakemake_profiling("profiling", project_name, args.type, args.output, env_path, args.profile, args.fraction, args.ani)

    ###
    # Dereplicating
    ###

    if args.command == "dereplicating":
        print(f"", flush=True)
        print(f"{HEADER1}STARTING DEREPLICATING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        bins_dir = args.bins_dir
        bins_file = args.bins_file

        if bins_dir and bins_file:
            print(f"Both bin path file and bin directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the bin info file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"DRAKKAR will run with the bins in the bin directory.")
            path_bins_to_json(bins_dir,args.output)
        else:
            print(f"No bin information was provided. DRAKKAR will try to guess the location of the bins.")
            if os.path.exists(f"{args.output}/cataloging/final/all_bin_paths.txt"):
                file_bins_to_json(f"{args.output}/cataloging/final/all_bin_paths.txt",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the cataloging module was run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an bin file (-B) or directory (-b).")
                return

        run_snakemake_dereplicating("dereplicating", project_name, args.output, env_path, args.profile, args.ani)

    ###
    # Annotating
    ###

    if args.command in ("annotating", "complete"):
        print(f"", flush=True)
        print(f"{HEADER1}STARTING ANNOTATING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        if args.command in ("annotating"):
            bins_dir = args.bins_dir
            bins_file = args.bins_file

        if args.command in ("complete"):
            bins_dir = None
            bins_file = None

        # Prepare bin dictionaries
        if bins_dir and bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(bins_dir,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        run_snakemake_annotating("annotating", project_name, args.annotation_type, args.output, env_path, args.profile)

    ###
    # Inspecting
    ###

    if args.command == "inspecting":
        print(f"", flush=True)
        print(f"{HEADER1}STARTING INSPECTING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(args.bins_dir,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(args.bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return
            
        run_snakemake_inspecting("annotating", project_name,  args.type, args.output, env_path, args.profile)


    ###
    # Expressing
    ###

    if args.command == "expressing":
        print(f"", flush=True)
        print(f"{HEADER1}STARTING EXPRESSING PIPELINE...{RESET}", flush=True)
        print(f"", flush=True)

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(args.bins_dir,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(args.bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        # Prepare read dictionaries

        reads_dir = args.reads_dir
        reads_file = args.reads_file

        if reads_dir and reads_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_transcriptome_to_json(reads_file,args.output)
        elif reads_file and not reads_dir:
            print(f"")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_transcriptome_to_json(reads_file,args.output)
        elif reads_dir and not reads_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reads from the input directory.")
            argument_transcriptome_to_json(reads_dir,args.output)
        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the reads.")
            if any(os.scandir(f"{args.output}/preprocessing/final")):
                argument_transcriptome_to_json(f"{args.output}/preprocessing/final",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        run_snakemake_inspecting("expressing", project_name, args.output, env_path, args.profile)


if __name__ == "__main__":
    main()
