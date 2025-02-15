import argparse
import os
import sys
import subprocess
import yaml
import pandas as pd
from pathlib import Path

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

def run_snakemake_complete(workflow, input_dir, output_dir):

    """ Run the complete workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} reads_dir={input_dir} output_dir={output_dir} reference={reference} cataloging_mode={mode}"
        f"--quiet rules"
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_preprocessing(workflow, input_dir, output_dir, reference, profile):

    """ Run the preprocessing workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} reads_dir={input_dir} output_dir={output_dir} reference={reference} "
        f"--quiet rules"
    ]

    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_cataloging(workflow, input_dir, output_dir, mode, profile):

    """ Run the cataloging workflow """

    snakemake_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config workflow={workflow} preprocess_dir={input_dir} output_dir={output_dir} cataloging_mode={mode} "
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


def display_drakkar():
    ascii_ship = r"""


                                    =+:   :##=   .###    +##    *#-
                                    -##-  -##*    ###-   +##    ##
                                      #*  .###    ###+   +##    ##+
                                      +#.  ###    ####   -##.   ###
                                      *#-  *##    *###   +##.   ##+
                                      *#=  =##    *###   -##-  .##+
                                      ##+  *##    ####   :##=  :##*              .-*##-
             *=#                      ##= .###    ####   =##+  .##*              +%%#--
               ++                    -##. .###    ####   +##*  .###              .%#
               #=                    ##*  :###    ####   +##*  .###               ++
              -%-                   ###.  +###    ####   =###   ###               =+:
              %#.                  ###+  .*##*    ###*   +###   ###               :++
             +@*                  ####   .###=   :###-   +###                      *+-
             @@*                 +###.   -###:   -###=   .##=                      @+=:
            *@@#:                        :+*#    =+=-                             =#++=
            @@@@%*=.                            .-                             :#@++*##-
           :@@@@@@++*+++*+*++-:.-*+..  =+--. .-.:-  ..  . .... .:=-::.+*%@#......+*****=
           :@@@@@@@#@%++#%%+==-===+---+####-:-+#*:+..:=*+*+:...:..:....:-==:..**++=+++=:
            #%%@**@@@@@@%#@@@#%@%%@-:-*%%##-..*#%##...=###*-..-+#===.+=+++---====*=--==.
            :%+*###%%%#*%%%%%%%#*%%%#%%%#*#######*#######*********++++*++++====-:::++=:
            :-*%@%%%#+#%%%%####*##%%%###########**#######**#####**=****+++++----====-:=
         -=.     =**#%%%%%%%%%*#%@@@@@@%##%%%%%%*#%%%%@@#**********+++=====++--:.       =
                 =-          =:    %@@@#*#######+*####@@:                    .
                            :=     %%--==.......:....:#@:
                                 :+@@--==:::::--.:=+*%@@:
    """

    ascii_text = r"""
     ██████████   ███████████     █████████   █████   ████ █████   ████   █████████   ███████████
    ░░███░░░░███ ░░███░░░░░███   ███░░░░░███ ░░███   ███░ ░░███   ███░   ███░░░░░███ ░░███░░░░░███
     ░███   ░░███ ░███    ░███  ░███    ░███  ░███  ███    ░███  ███    ░███    ░███  ░███    ░███
     ░███    ░███ ░██████████   ░███████████  ░███████     ░███████     ░███████████  ░██████████
     ░███    ░███ ░███░░░░░███  ░███░░░░░███  ░███░░███    ░███░░███    ░███░░░░░███  ░███░░░░░███
     ░███    ███  ░███    ░███  ░███    ░███  ░███ ░░███   ░███ ░░███   ░███    ░███  ░███    ░███
     ██████████   █████   █████ █████   █████ █████ ░░████ █████ ░░████ █████   █████ █████   █████
    ░░░░░░░░░░   ░░░░░   ░░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░ ░░░░░   ░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░░
    """

    ascii_intro = r"""

    By Antton Alberdi [antton.alberdi@sund.ku.dk]
    DRAKKAR is a snakemake-based genome-resolved metagenomics pipeline optimised for Mjolnir.
    """

    print(ascii_ship)
    print(ascii_text)
    print(ascii_intro)

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

    ###
    # Check screen session
    ###

    if "STY" not in os.environ:
        print(" ⚠️   WARNING: You are not running this script inside a 'screen' session.")
        print("     Running long processes outside of 'screen' may cause issues if your session is disconnected.")
        print("     \n📌 To start a screen session, use:  screen -S mysession")
        print("         Then run this script inside the screen session.\n")

        # Prompt user to continue
        while True:
            user_input = input("    👉 Type '1' to ignore this warning and continue, or Ctrl+C to exit: ")
            if user_input.strip() == "1":
                break  # Continue execution
            else:
                print("     Invalid input. Please type '1' to continue.")

    ###
    # Processing of sample detail file
    ###

    if args.file:
        file_path = Path(args.file).resolve()
        if not file_path.exists():
            print(f"ERROR: Sample info file '{file_path}' was not found.")
            return
        try:
            df = pd.read_csv(file_path, sep="\t")
            if "sample" not in df.columns:
                print(f"WARNING: The mandatory column 'sample' was not found in the file.")
                unique_samples = "N/A"
            else:

                ####
                # Process samples
                ####

                unique_samples = df["sample"].nunique()

                FILE_TO_SAMPLE = dict(zip(df["sample"], df["reference"]))
                print(FILE_TO_SAMPLE)

                total_datafiles = len(df)
                print(f"")
                print(f" Running DRAKKAR with {total_datafiles} files belonging to {unique_samples} samples.")
                print(f" Preparing input data based on the sample info file.")

                INFOFILE = True

                # Ensure required columns exist
                required_columns = {"sample", "rawreads1", "rawreads2"}
                if not required_columns.issubset(df.columns):
                    raise ValueError(f"Missing required columns in {TABLE_FILE}. Expected: {required_columns}")

                # Group file paths by sample
                samples = {}
                for _, row in df.iterrows():
                    sample = row["sample"]
                    rawreads1, rawreads2 = row["rawreads1"], row["rawreads2"]

                    if sample not in samples:
                        samples[sample] = {"rawreads1": [], "rawreads2": []}

                    samples[sample]["rawreads1"].append(rawreads1)
                    samples[sample]["rawreads2"].append(rawreads2)

                # Create output directories and concatenate files
                for sample, reads in samples.items():
                    sample_dir = Path(args.output).resolve() / "data" / "reads"
                    sample_dir.mkdir(parents=True, exist_ok=True)  # Create output directory

                    rawreads1_output = sample_dir / f"{sample}_1.fq.gz"
                    rawreads2_output = sample_dir / f"{sample}_2.fq.gz"

                    print(f"   📂 Processing sample: {sample}")

                    # Concatenate rawreads1 using `cat`
                    rawreads1_cmd = f"cat {' '.join(reads['rawreads1'])} > {rawreads1_output}"
                    rawreads2_cmd = f"cat {' '.join(reads['rawreads2'])} > {rawreads2_output}"

                    subprocess.run(rawreads1_cmd, shell=True, check=True)
                    subprocess.run(rawreads2_cmd, shell=True, check=True)

                    INPUTPATH=sample_dir


                ####
                # Process reference genomes
                ####

                # Extract mapping of samples to references
                SAMPLE_TO_REFERENCE = dict(zip(df["sample"], df["reference"]))
                print(SAMPLE_TO_REFERENCE)

        except Exception as e:
            print(f"Error reading file: {e}")

    else:
        if args.input:
            print(f"")
            print(f"    No sample info file was provided. Drakkar will guess samples from the provided directory.")
            INPUTPATH=Path(args.input).resolve()
        else:
            print(f"")
            print(f"Please provide either an input directory (-i) or a sample info file (-f)")
            return  # Exit after file processing

    ###
    # Launch snakemake commands
    ###
    print(f"")
    print(f"Starting Snakemake pipeline...")
    print(f"")

    # Relative paths are turned into absolute paths
    if args.command == "complete":
        run_snakemake_complete(args.command, args.input, args.output, args.reference)
    elif args.command == "preprocessing":
        run_snakemake_preprocessing(args.command, INPUTPATH, Path(args.output).resolve(), Path(args.reference).resolve() if args.reference else None, args.profile if args.profile else "slurm")
    elif args.command == "cataloging":
        run_snakemake_cataloging(args.command, Path(args.input).resolve(), Path(args.output).resolve(), args.mode if args.mode else ["individual"], args.profile if args.profile else "slurm")
    elif args.command == "annotation":
        run_snakemake_annotation(args.command, args.assembly, args.output)
    elif args.command == "quantification":
        run_snakemake_quantification(args.command, args.assembly, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
