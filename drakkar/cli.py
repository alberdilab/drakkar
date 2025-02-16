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

def run_snakemake_complete(workflow, input_dir, output_dir, reference, mode, profile):

    """ Run the complete workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} output_dir={output_dir} reference={reference} cataloging_mode={mode} "
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_preprocessing(workflow, output_dir, reference, profile):

    """ Run the preprocessing workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} output_dir={output_dir} reference={reference} "
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_cataloging(workflow, output_dir, mode, profile):

    """ Run the cataloging workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} output_dir={output_dir} cataloging_mode={mode} "
        f"--quiet rules"
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_dereplicating(workflow, assembly_dir, output_dir):
    """ Run the dereplicating workflow """
    snakemake_command = [
        f"module load {SNAKEMAKE_MODULE}",
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config",
            f"workflow={workflow}",
            f"assembly_dir={assembly_dir}",
            f"output_dir={output_dir}"
    ]
    subprocess.run(" && ".join(snakemake_command), shell=True, check=True)

def run_snakemake_annotation(workflow, assembly_dir, output_dir):
    """ Run the annotation workflow """
    snakemake_command = [
        f"module load {SNAKEMAKE_MODULE}",
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config",
            f"workflow={workflow}",
            f"assembly_dir={assembly_dir}",
            f"output_dir={output_dir}"
    ]
    subprocess.run(" && ".join(snakemake_command), shell=True, check=True)

def run_snakemake_quantification(workflow, assembly_dir, output_dir):
    """ Run the quantification workflow """
    snakemake_command = [
        f"module load {SNAKEMAKE_MODULE}",
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config",
            f"workflow={workflow}",
            f"assembly_dir={assembly_dir}",
            f"output_dir={output_dir}"
    ]
    subprocess.run(" && ".join(snakemake_command), shell=True, check=True)

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
    subparser_complete.add_argument("-o", "--output", required=True, help="Output directory")
    subparser_complete.add_argument("-r", "--reference", required=False, help="Reference host genome")

    subparser_preprocessing = subparsers.add_parser("preprocessing", help="Run the preprocessing workflow (quality-filtering and host removal)")
    subparser_preprocessing.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_preprocessing.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_preprocessing.add_argument("-o", "--output", required=True, help="Output directory")
    subparser_preprocessing.add_argument("-r", "--reference", required=False, help="Reference host genome file (fna.gz)")
    subparser_preprocessing.add_argument("-p", "--profile", required=False, help="Snakemake profile")

    subparser_cataloging = subparsers.add_parser("cataloging", help="Run the cataloging workflow (assemly and binning)")
    subparser_cataloging.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_cataloging.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_cataloging.add_argument("-o", "--output", required=True, help="Output directory")
    subparser_cataloging.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,custom,all)")
    subparser_cataloging.add_argument("-p", "--profile", required=False, help="Snakemake profile")

    subparser_dereplication = subparsers.add_parser("dereplication", help="Run the dereplication workflow")
    subparser_dereplication.add_argument("-b", "--bins", required=True, help="Bins directory")
    subparser_dereplication.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_annotation = subparsers.add_parser("annotation", help="Run the annotation workflow")
    subparser_annotation.add_argument("-a", "--assembly", required=True, help="Assembly directory")
    subparser_annotation.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_quantification = subparsers.add_parser("quantification", help="Run the quantification workflow")
    subparser_quantification.add_argument("-a", "--assembly", required=True, help="Assembly directory")
    subparser_quantification.add_argument("-o", "--output", required=True, help="Output directory")

    args = parser.parse_args()

    # Display ASCII logo before running any command or showing help
    display_drakkar()

    # Check screen session
    check_screen_session()

    ###
    # Preprocessing
    ###

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

        elif args.file and not args.reference:
            if check_reference_columns(args.file):
                print(f"")
                print(f"DRAKKAR will extract the reference genome information from the sample info file.")
                file_references_to_json(args.file,args.output)
                REFERENCE = True

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

            with open(f"{args.output}/data/preprocessed_to_reads1.json", "r") as f:
                PREPROCESSED_TO_READS1 = json.load(f)

        elif args.file and not args.input:
            print(f"")
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

            with open(f"{args.output}/data/sample_to_reads1.json", "r") as f:
                PREPROCESSED_TO_READS1 = json.load(f)

        elif args.input and not args.file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)

            with open(f"{args.output}/data/preprocessed_to_reads1.json", "r") as f:
                PREPROCESSED_TO_READS1 = json.load(f)

        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the preprocessing data")
            if any(os.scandir(f"{args.output}/preprocessed/final")):
                argument_preprocessed_to_json(f"{args.output}/preprocessed/final",args.output)

                with open(f"{args.output}/data/preprocessed_to_reads1.json", "r") as f:
                    PREPROCESSED_TO_READS1 = json.load(f)
            else:
                print(f"ERROR: No preprocessed data was found in the output directory.")
                print(f"    Please, provide an input directory or sample info file to proceed")
                return

        # Generate the assembly dictionary

        samples = list(PREPROCESSED_TO_READS1.keys())

        if args.mode == "individual":
            INDIVIDUAL_MODE=True
            ALL_MODE=False
        elif args.mode == "all":
            INDIVIDUAL_MODE=False
            ALL_MODE=True
        else:
            INDIVIDUAL_MODE=False
            ALL_MODE=True

        file_assemblies_to_json(args.input,samples,INDIVIDUAL_MODE,ALL_MODE,args.output)

    ###
    # Launch snakemake commands
    ###
    print(f"")
    print(f"Starting Snakemake pipeline...")
    print(f"")

    # Relative paths are turned into absolute paths
    if args.command == "complete":
        run_snakemake_complete(args.command, args.input, args.output, args.reference, args.mode if args.mode else ["individual"], args.profile if args.profile else "slurm")
    elif args.command == "preprocessing":
        run_snakemake_preprocessing(args.command, Path(args.output).resolve(), REFERENCE, args.profile if args.profile else "slurm")
        display_end()
    elif args.command == "cataloging":
        run_snakemake_cataloging(args.command, Path(args.output).resolve(), args.mode if args.mode else ["individual"], args.profile if args.profile else "slurm")
    elif args.command == "annotation":
        run_snakemake_annotation(args.command, args.assembly, args.output)
    elif args.command == "quantification":
        run_snakemake_quantification(args.command, args.assembly, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
