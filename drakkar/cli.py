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
# Define workflow launching functions
###

def run_unlock(workflow, output_dir, profile):

    unlock_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)

def run_snakemake_complete(workflow, project_name, input_dir, output_dir, reference, mode, profile):

    """ Run the complete workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} reference={reference} "
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_preprocessing(workflow, project_name, output_dir, reference, profile):

    """ Run the preprocessing workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} reference={reference} "
    ]

    try:
        subprocess.run(snakemake_command, shell=False, check=True)
    except subprocess.CalledProcessError as e:
        error_message = str(e)
        if "Error: Directory cannot be locked" in error_message:
            display_unlock()
        else:
            print(f"\nERROR: Snakemake failed with exit code {e.returncode}!", file=sys.stderr)
            print(f"ERROR: Check the Snakemake logs for more details.", file=sys.stderr)
            sys.exit(1)


def run_snakemake_cataloging(workflow, project_name, output_dir, profile):

    """ Run the cataloging workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_profiling(workflow, project_name, profiling_type, output_dir, profile):
    """ Run the profiling workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} profiling_type={profiling_type} output_dir={output_dir}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_annotating(workflow, project_name, output_dir, profile):
    """ Run the profiling workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir}"
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
    subparser_complete.add_argument("-i", "--input", required=True, help="Input directory")
    subparser_complete.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_complete.add_argument("-r", "--reference", required=False, help="Reference host genome")
    subparser_complete.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_complete.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    subparser_preprocessing = subparsers.add_parser("preprocessing", help="Run the preprocessing workflow (quality-filtering and host removal)")
    subparser_preprocessing.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_preprocessing.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_preprocessing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_preprocessing.add_argument("-r", "--reference", required=False, help="Reference host genome file (fna.gz)")
    subparser_preprocessing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    subparser_cataloging = subparsers.add_parser("cataloging", help="Run the cataloging workflow (assemly and binning)")
    subparser_cataloging.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_cataloging.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_cataloging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_cataloging.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_cataloging.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    subparser_profiling = subparsers.add_parser("profiling", help="Run the profiling workflow")
    subparser_profiling.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_profiling.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_profiling.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_profiling.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_profiling.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_profiling.add_argument("-t", "--type", required=False, default="genomes", help="Either genomes or pangenomes profiling type. Default: genomes")
    subparser_profiling.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    subparser_annotating = subparsers.add_parser("annotating", help="Run the cataloging workflow (assemly and binning)")
    subparser_annotating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_annotating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_annotating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_annotating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    subparser_unlock = subparsers.add_parser("unlock", help="Unlock snakemake")
    subparser_unlock.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_unlock.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile")

    args = parser.parse_args()

    # Display ASCII logo before running any command or showing help
    display_drakkar()

    # Check screen session
    check_screen_session()

    # Extract project name from output path


    ###
    # Preprocessing
    ###

    if args.command == "unlock":
        run_unlock(args.command, args.output, args.profile)
    else:
        project_name = os.path.basename(os.path.normpath(args.output))

    if args.command == "preprocessing":

        # Generate raw data dictionaries

        if args.file and args.input:
            print(f"")
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.file and not args.input:
            print(f"")
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.input and not args.file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_samples_to_json(args.input,args.output)

        else:
            print(f"")
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
            if check_reference_columns(args.file):
                print(f"")
                print(f"No sample info file was provided.")
                print(f"DRAKKAR will use the reference genome file.")
                argument_references_to_json(args.reference,f"{args.output}/data/sample_to_reads1.json",args.output)
                REFERENCE = True

        else:
            print(f"")
            print(f"Running DRAKKAR without mapping against a reference genome")
            REFERENCE = False

    ###
    # Cataloging
    ###

    if args.command == "cataloging":

        # Generate preprocessing data dictionaries

        if args.file and args.input:
            print(f"")
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        elif args.file and not args.input:
            print(f"")
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_preprocessed_to_json(args.file,args.output)
        elif args.input and not args.file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        else:
            print(f"")
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


    ###
    # Profiling
    ###

    if args.command == "profiling":

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_bins_to_json(args.input,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"")
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_bins_to_json(args.file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_bins_to_json(args.input,args.output)
        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/cataloging/final/all_bin_paths.txt"):
                file_bins_to_json(f"{args.output}/cataloging/final/all_bin_paths.txt",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        # Prepare read dictionaries
        if args.reads_dir and args.reads_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_bins_to_json(args.input,args.output)
        elif args.reads_file and not args.reads_dir:
            print(f"")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_bins_to_json(args.file,args.output)
        elif args.reads_dir and not args.reads_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reads from the input directory.")
            path_bins_to_json(args.input,args.output)
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

    ###
    # Annotating
    ###

    if args.command == "annotating":

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(args.input,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"")
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(args.file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(args.input,args.output)
        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

    ###
    # Launch snakemake commands
    ###
    print(f"")
    print(f"Starting Snakemake pipeline...")
    print(f"")

    # Relative paths are turned into absolute paths
    if args.command == "complete":
        run_snakemake_complete(args.command, project_name, args.input, args.output, args.reference, args.profile)
        display_end()
    elif args.command == "preprocessing":
        run_snakemake_preprocessing(args.command, project_name, Path(args.output).resolve(), REFERENCE, args.profile)
        display_end()
        preprocessing_summary(f"{args.output}/preprocessing.tsv")
    elif args.command == "cataloging":
        run_snakemake_cataloging(args.command, project_name, Path(args.output).resolve(), args.profile)
        display_end()
    elif args.command == "profiling":
        run_snakemake_profiling(args.command, project_name, args.type, args.output, args.profile)
        display_end()
    elif args.command == "annotating":
        run_snakemake_annotating(args.command, project_name, args.output, args.profile)
        display_end()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
