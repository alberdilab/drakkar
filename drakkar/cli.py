import argparse
import os
import subprocess
from pathlib import Path
from drakkar.config import SNAKEMAKE_MODULE  # Load from config.py

###
# Define workflow launching functions
###

def run_snakemake_complete(workflow, input_dir, output_dir):
    """ Run the complete workflow """
    snakemake_command = [
        f"module load {SNAKEMAKE_MODULE}",
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config",
            f"workflow={workflow}",
            f"reads_dir={input_dir}",
            f"output_dir={output_dir}"
    ]
    subprocess.run(" && ".join(snakemake_command), shell=True, check=True)

def run_snakemake_preprocessing(workflow, input_dir, output_dir):
    """ Run the preprocessing workflow """
    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {SNAKEMAKE_MODULE} && "
        "snakemake "
        f"-s {str(Path(__file__).parent / 'workflow' / 'Snakefile')} "
        f"--config workflow={workflow} reads_dir={input_dir} output_dir={output_dir}"
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_assembly(workflow, input_dir, output_dir):
    """ Run the assembly workflow """
    snakemake_command = [
        f"module load {SNAKEMAKE_MODULE}",
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config",
            f"workflow={workflow}",
            f"reads_dir={input_dir}",
            f"output_dir={output_dir}"
    ]
    subprocess.run(" && ".join(snakemake_command), shell=True, check=True)

def run_snakemake_binning(workflow, assembly_dir, output_dir):
    """ Run the binning workflow """
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
    parser = argparse.ArgumentParser(description="Drakkar: A Snakemake-based workflow for sequencing analysis")
    subparsers = parser.add_subparsers(dest="command")

    # Define subcommands for each workflow
    subparser_complete = subparsers.add_parser("complete", help="Run the complete workflow")
    subparser_complete.add_argument("-i", "--input", required=True, help="Input directory")
    subparser_complete.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_assembly = subparsers.add_parser("preprocessing", help="Run the preprocessing workflow")
    subparser_assembly.add_argument("-i", "--input", required=True, help="Input directory")
    subparser_assembly.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_assembly = subparsers.add_parser("assembly", help="Run the assembly workflow")
    subparser_assembly.add_argument("-i", "--input", required=True, help="Input directory")
    subparser_assembly.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_binning = subparsers.add_parser("binning", help="Run the binning workflow")
    subparser_binning.add_argument("-a", "--assembly", required=True, help="Assembly directory")
    subparser_binning.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_annotation = subparsers.add_parser("annotation", help="Run the annotation workflow")
    subparser_annotation.add_argument("-a", "--assembly", required=True, help="Assembly directory")
    subparser_annotation.add_argument("-o", "--output", required=True, help="Output directory")

    subparser_quantification = subparsers.add_parser("quantification", help="Run the quantification workflow")
    subparser_quantification.add_argument("-a", "--assembly", required=True, help="Assembly directory")
    subparser_quantification.add_argument("-o", "--output", required=True, help="Output directory")

    args = parser.parse_args()

    if args.command == "complete":
        run_snakemake_complete(args.command, args.input, args.output)
    elif args.command == "preprocessing":
        run_snakemake_preprocessing(args.command, args.input, args.output)
    elif args.command == "assembly":
        run_snakemake_assembly(args.command, args.input, args.output)
    elif args.command == "binning":
        run_snakemake_binning(args.command, args.assembly, args.output)
    elif args.command == "annotation":
        run_snakemake_annotation(args.command, args.assembly, args.output)
    elif args.command == "quantification":
        run_snakemake_quantification(args.command, args.assembly, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
