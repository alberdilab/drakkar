import subprocess
import sys
from pathlib import Path

from drakkar.cli_context import CONFIG_PATH, DEFAULT_CATALOGING_BINNERS, PACKAGE_DIR, config_vars
from drakkar.cli_validation import default_resource_args, resource_config
from drakkar.display import display_end, display_unlock
from drakkar.output import print
from drakkar.run_metadata import run_subprocess_with_logging

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

def run_snakemake_environments(workflow, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    cmd = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {Path.cwd()} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} workflow={workflow} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(cmd, run_info=run_info, workflow_name=workflow)

def run_snakemake_preprocessing(
    workflow,
    project_name,
    output_dir,
    reference,
    env_path,
    profile,
    fraction=False,
    nonpareil=False,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):

    """ Run the preprocessing workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} reference={reference} fraction={fraction} nonpareil={nonpareil} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
        f"--slurm-delete-logfiles-older-than 0"
    ]

    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_cataloging(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
    binners=DEFAULT_CATALOGING_BINNERS,
):

    """ Run the cataloging workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} binners={binners} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
    ]

    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_cataloging2(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    memory_multiplier=1,
    time_multiplier=1,
    binners=DEFAULT_CATALOGING_BINNERS,
):

    """ Run the cataloging workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} binners={binners} {resource_overrides}"
        f"{default_resources}"
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

def run_snakemake_profiling(
    workflow,
    project_name,
    profiling_type,
    output_dir,
    env_path,
    profile,
    fraction,
    ani,
    ignore_quality,
    quality_file,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the profiling workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} profiling_type={profiling_type} output_dir={output_dir} fraction={fraction} DREP_ANI={ani} "
        f"IGNORE_QUALITY={ignore_quality} QUALITY_FILE={bool(quality_file)} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_dereplicating(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    ani,
    ignore_quality,
    quality_file,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the dereplicating workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} DREP_ANI={ani} IGNORE_QUALITY={ignore_quality} QUALITY_FILE={bool(quality_file)} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_annotating(
    workflow,
    project_name,
    annotating_type,
    output_dir,
    env_path,
    profile,
    gtdb_version=None,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the profiling workflow """

    gtdb_config = f"gtdb_version={gtdb_version} " if gtdb_version else ""
    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} annotating_type={annotating_type} output_dir={output_dir} {gtdb_config}{resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_inspecting(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    """ Run the profiling workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c", 
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_expressing(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    """ Run the expressing workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_database(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    database_name,
    database_directory,
    database_version,
    download_runtime,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """Run a single database preparation workflow."""

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"database_name={database_name} database_directory={database_directory} database_version={database_version} "
        f"database_download_runtime={download_runtime} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)
