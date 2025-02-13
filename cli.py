import argparse
import os
import subprocess
from pathlib import Path

def run_snakemake(workflow, input_dir, output_dir, assembly_dir=None):
    snakemake_command = [
        "snakemake",
        "-s", str(Path(__file__).parent / "workflow" / "Snakefile"),
        "--config", f"input_dir={input_dir}", f"output_dir={output_dir}"
    ]
    if assembly_dir:
        snakemake_command.append(f"assembly_dir={assembly_dir}")

    if workflow != "complete":
        snakemake_command.extend(["--rulegraph", workflow])

    subprocess.run(snakemake_command, check=True)

def main():
    parser = argparse.ArgumentParser(description="Drakkar: A Snakemake-based workflow for sequencing analysis")
    subparsers = parser.add_subparsers(dest="command")

    workflows = ["complete", "assembly", "binning", "annotation", "quantification"]
    for workflow in workflows:
        subparser = subparsers.add_parser(workflow, help=f"Run {workflow} workflow")
        subparser.add_argument("-i", "--input", required=(workflow in ["assembly", "complete"]), help="Input directory")
        subparser.add_argument("-o", "--output", required=True, help="Output directory")
        if workflow in ["binning", "annotation", "quantification"]:
            subparser.add_argument("-a", "--assembly", required=True, help="Assembly directory")

    args = parser.parse_args()

    if args.command:
        run_snakemake(args.command, args.input, args.output, getattr(args, 'assembly', None))
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
