import os
import shutil
import sys
from pathlib import Path

from drakkar.cli_context import ERROR, INFO, RESET
from drakkar.downloads import is_url
from drakkar.output import print, prompt
from drakkar.system_checks import is_snakemake_locked

def overwrite_output_directory(output_dir):
    output_path = Path(output_dir)
    if output_path.exists():
        shutil.rmtree(output_path)

def prompt_overwrite_locked_directory(output_dir):
    print(f"{INFO}INFO:{RESET} The output directory is locked: {output_dir}")
    print("This usually means a previous Snakemake run ended unexpectedly.")
    print(f"Use 'drakkar logging -o {output_dir}' to inspect the latest Snakemake log before unlocking or overwriting.")
    print("Overwriting will delete the entire directory and rerun from scratch.")
    while True:
        response = prompt("Delete the locked directory and overwrite it? [y/N]: ").strip().lower()
        if response in {"y", "yes"}:
            return True
        if response in {"", "n", "no"}:
            return False
        print("Please answer y or n.")

def prepare_output_directory(output_dir, overwrite=False):
    output_path = Path(output_dir)
    if not is_snakemake_locked(str(output_path)):
        return True

    if overwrite:
        print(f"{INFO}INFO:{RESET} Removing locked output directory and overwriting: {output_path}")
        overwrite_output_directory(output_path)
        return True

    if sys.stdin.isatty():
        if prompt_overwrite_locked_directory(output_path):
            overwrite_output_directory(output_path)
            return True
    else:
        print(f"{ERROR}ERROR:{RESET} Output directory is locked and no interactive prompt is available.")

    print(f"{ERROR}ERROR:{RESET} Output directory remains locked: {output_path}")
    print(f"{INFO}Use 'drakkar logging -o {output_path}' to inspect the latest Snakemake log, --overwrite to delete it automatically, or 'drakkar unlock -o {output_path}' if you only want to clear the Snakemake lock.{RESET}")
    return False

def validate_path(path_value, label, expect_dir=False, allow_url=False):
    if not path_value:
        return True
    if allow_url and not expect_dir and is_url(path_value):
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

def validate_launch_metadata_directory(output_dir):
    output_path = Path(output_dir)
    existing_path = output_path
    while not existing_path.exists() and existing_path.parent != existing_path:
        existing_path = existing_path.parent

    if not existing_path.exists():
        print(f"{ERROR}ERROR:{RESET} Cannot find a parent directory for output path: {output_path}")
        return False
    if not existing_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path parent is not a directory: {existing_path}")
        return False

    access_target = output_path if output_path.exists() else existing_path
    if output_path.exists() and not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path is not a directory: {output_path}")
        return False
    if not os.access(access_target, os.R_OK | os.W_OK | os.X_OK):
        print(f"{ERROR}ERROR:{RESET} Cannot write Drakkar run metadata in: {output_path}")
        print("The current/output directory must be readable, writable, and searchable.")
        print("Run drakkar from a writable directory or pass -o/--output to a writable output directory.")
        return False

    return True
