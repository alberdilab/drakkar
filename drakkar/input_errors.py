import os
import sys

from drakkar.output import print

class InputFileError(Exception):
    pass

class DownloadError(InputFileError):
    pass

def require_non_empty_file(file_path, description):
    if not os.path.isfile(file_path):
        raise InputFileError(f"{description} not found: {file_path}")
    if os.path.getsize(file_path) == 0:
        raise InputFileError(f"{description} is empty: {file_path}")
    return file_path

def report_input_resolution_errors(errors):
    print("ERROR: DRAKKAR stopped before launching Snakemake because input files could not be prepared.")
    print("Problems detected:")
    for error in errors:
        print(f"  - {error}")
    print("Fix the listed paths or URLs and rerun the same command. Snakemake was not launched.")
    sys.exit(1)
