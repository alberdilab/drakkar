import os
import sys
import yaml
import re
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict

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
    DRAKKAR is a snakemake-based genome-resolved metagenomics pipeline optimised for Mjolnir. Snakemake
    works along with Slurm to conduct the long pipeline using the optimal memory and time resources.
    This means that while running DRAKKAR, many Slurm jobs will be automatically sent to the Mjolnir queue.
    Depending on how busy the server is, the pipeline will run faster or slower. While running, Snakemake will
    show the progress on screen. DRAKKAR has multiple usage modes, so make sure you read the documentation
    properly before using it in https://github.com/alberdilab/drakkar.
    """

    print(ascii_ship)
    print(ascii_text)
    print(ascii_intro)

def display_end():
    ascii_dragon = r"""

                                         =
                                       .=*.-+
                                  .:----+**+==-.
                         :-: .:-===++====*#*==+*+:
                         :*#*****=+++++*=+++=+++++
                             .=***+++=----=+**#*=--=.
                              -+++=--=+**++=++=++===+=+++-
                             -+++=-=****+***=+*++=#**++**+
                            .=++=:+##*+*+=+++*#%#-=++***-
                            :=+-:=####*+===+*##+++   +#.
                           .+#+--+#####+-==++## :+.   .
                           :-=+*-+####     =+**:   -
                           :=++++=+**-      :+**+:
                           -+***+=+**+.       -:
                           =++**+****+
                           =+++***##=
                           =*++***##=             CONGRATULATIONS!
                           =*+*+*###+             DRAKKAR finished successfully.
                           =*++***##+             Check the output directory to
                           =*+****#**             find your results.
                           -++*****#*
                           -+=**+###*.            If you don't find what you need
                           :++++***##:            read the documentation in
                           .+++******+            https://github.com/alberdilab/drakkar
                            =+*******+.
                            -*+******+-
                            :+*+****+++.
                             =+++*++=++=

    """
    print(ascii_dragon)

def check_screen_session():
    """Checks if the script is running inside a 'screen' session. If not, warns the user and asks for confirmation."""
    if "STY" not in os.environ:
        print("\n ⚠️   WARNING: You are not running this script inside a 'screen' session.")
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

def check_reference_columns(file_path):
    """Checks if a file contains 'reference_name' and 'reference_path' columns with non-NA values."""
    # Read the file (assumed to be TSV, change sep="," for CSV)
    df = pd.read_csv(file_path, sep="\t")
    # Check if required columns exist
    required_columns = {"reference_name", "reference_path"}
    if not required_columns.issubset(df.columns):
        return False
    # Check if both columns have at least one non-NA value
    if df["reference_name"].dropna().empty or df["reference_path"].dropna().empty:
        return False
    return True

def check_assembly_column(file_path):
    """Checks if a file contains the column 'assembly' with non-NA values."""
    # Read the file (assumed to be TSV, change sep="," for CSV)
    df = pd.read_csv(file_path, sep="\t")
    # Check if required columns exist
    required_columns = {"assembly"}
    if not required_columns.issubset(df.columns):
        return False
    # Check if both columns have at least one non-NA value
    if df["assembly"].dropna().empty:
        return False
    return True

def file_samples_to_json(infofile, output):
    # Load sample info file
    df = pd.read_csv(infofile, sep="\t")

    # Initialize dictionaries with lists
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Populate the dictionaries
    for _, row in df.iterrows():
        SAMPLE_TO_READS1[row["sample"]].append(str(Path(row["rawreads1"]).resolve()))
        SAMPLE_TO_READS2[row["sample"]].append(str(Path(row["rawreads2"]).resolve()))

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    # Save dictionaries to JSON files
    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)


def file_references_to_json(infofile, output):
    # Load sample info file
    df = pd.read_csv(infofile, sep="\t")

    REFERENCE_TO_FILE = {ref_name: str(Path(ref_path).resolve()) for ref_name, ref_path in zip(df["reference_name"], df["reference_path"])}
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    SAMPLE_TO_REFERENCE = dict(zip(df["sample"], df["reference_name"]))
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

#def file_assemblies_to_json(infofile):

def argument_samples_to_json(argument, output):
    # Define the directory containing the raw reads
    READS_DIR = Path(argument).resolve()

    # Initialize dictionaries
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(READS_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(READS_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    SAMPLE_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    SAMPLE_TO_READS2[sample_name].append(full_path)

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_references_to_json(argument, sample_to_reads, output):
    REFERENCE_TO_FILE = {"reference": [f"{argument}"]}
    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    with open(f"{sample_to_reads}", "r") as f:
        SAMPLE_TO_READS = json.load(f)

    SAMPLE_TO_REFERENCE = {sample: "reference" for sample in SAMPLE_TO_READS.keys()}
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

def file_preprocessed_to_json(infofile, output):
    # Load sample info file
    df = pd.read_csv(infofile, sep="\t")

    # Initialize dictionaries with lists
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Populate the dictionaries
    for _, row in df.iterrows():
        SAMPLE_TO_READS1[row["sample"]].append(str(Path(row["rawreads1"]).resolve()))
        SAMPLE_TO_READS2[row["sample"]].append(str(Path(row["rawreads2"]).resolve()))

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    # Save dictionaries to JSON files
    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/preprocessed_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/preprocessed_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_preprocessed_to_json(argument, output):
    # Define the directory containing the raw reads
    PREPROCESSED_DIR = Path(argument).resolve()

    # Initialize dictionaries
    PREPROCESSED_TO_READS1 = defaultdict(list)
    PREPROCESSED_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(PREPROCESSED_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(PREPROCESSED_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    PREPROCESSED_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    PREPROCESSED_TO_READS2[sample_name].append(full_path)

    # Convert defaultdict to standard dict (optional)
    PREPROCESSED_TO_READS1 = dict(PREPROCESSED_TO_READS1)
    PREPROCESSED_TO_READS2 = dict(PREPROCESSED_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/preprocessed_to_reads1.json", "w") as f:
        json.dump(PREPROCESSED_TO_READS1, f)

    with open(f"{output}/data/preprocessed_to_reads2.json", "w") as f:
        json.dump(PREPROCESSED_TO_READS2, f)

def file_assemblies_to_json(infofile=None, samples=None, individual=False, all=False, output=False):

    assemblies = defaultdict(list)

    if infofile is not None:
        df = pd.read_csv(infofile, sep="\t")
        for _, row in df.iterrows():
            sample = row['sample']
            assembly_list = row['assembly'].split(',')

            for assembly in assembly_list:
                assemblies[assembly].append(sample)

    if samples:
        if individual:
            for sample in samples:
                assemblies[sample].append(sample)

        if all:
            assemblies["all"].extend(samples)

    ASSEMBLY_TO_SAMPLE = {key: list(set(value)) for key, value in assemblies.items()} # remove duplicates

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/assembly_to_samples.json", "w") as f:
        json.dump(ASSEMBLY_TO_SAMPLE, f)
