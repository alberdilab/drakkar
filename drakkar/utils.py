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
    Source code: https://github.com/alberdilab/drakkar
    Tutorial: https://drakkar.readthedocs.io/
    """

    print(ascii_ship)
    print(ascii_text)
    print(ascii_intro)

def display_unlock():
    ascii_swords = r"""
    @@%*:.                  DRAKKAR OUTPUT DIRECTORY IS LOCKED              .:*%@@
    =@@@%@@*:.      -------------------------------------------------     .-*@@%@@@=
    .#@#@*:+%@#-.         Probably because a DRAKKAR run ended          ..-%@%+-#@#@#.
     .@%:#@+..-%@#:.       unexpectedly. Use 'squeue' to check        .:%@%=::*@%-%@.
      :@#:-%@=..:=@@+.        whether any slurm job is still        .+@@+-::+@%--%@:
       :@#:.-@@-...-#@#:.      running. Wait until finishing      .:#@#-:::=@@=::%@:.
        :%%:..=@%-....+@@-.       or kill the job using         ..-@@+::::-%@+::-%%:
         .#@-...+@%:....-@@=       'scancel'. Then, use       .=@@=::::-%@+:::-@#.
          .*@+...:*@*:....-%@+.      'drakkar unlock'      .+@%=::::-*@*-:::+@*.
            -@%:...:#@+:....-%@+.      to unlock the      .+@%-:::::+@#-:::-%@-
             .%@=....:%@=.....:%@=.     directory.     ..=@%-:::::=@%-::::=@%.
              .=@#.....-@@:.....-%@-                 .-@%-:::::-@@=:::::#@=.
                .%@=....:=@#:.....-%@-             .-@%=:::::-%@+:::::=@%.
                  -@%:....:*@*:.....-@%:         .:%@=:::::-#@*:::::-%@-.
                   .*@*:....:#@+......=@%.     ..%@+::::::*@#-:::::*@*.
                    .:@@=.....:%@=.....:+@*.  .*@*::::::=@%-:::::=@@:.
                      .=@%:.....-%@-.....:#@++@#-:::::=@@=:::::-%@=.
                        .*@*:.....-@%-.....-%@+:::::-%@=:::::-#@*.
                          .#@+......=@%:.....+@@-:-%@+::::::+@#.
                            :%@=.....:+@*:....:#@%@*::::::=@%:.
                             .:%@-.....:#@+.....:%@+::::=@%-.
                               .-@%-.....:%@-.....-@%-=@@-.
                                  -@%:.....:@@:.....=@@-
                                .+@#*@%:.....=@%:....:#@+.
                              .-@@=:::*@%:.....+@*:....-%@-.
                        ..  ..%@*:::::*@%@%:....:#@+:....+@%..   ..
                     .:##:..*@#:::::+@%-::*@#:....:%@=....:#@*..:##:..
                    :%%++@@@%-::::=@%-::::-#@@*:....-%@-....:%@@%++%%:
                    -@*---=@@=::-%@=::::-#@#..#@*:....-@%-..-%@=---*@-
                    .=@%----+@%%@+:::::#@%.   ..%@*:....=@%%@+----%@=.
                      :%@#----*@#-:::*@%:.      .:%@*...:#@*----#@%.
                    .*@%#%@*----#@*+@%:.          .:%@+*@#----*@%#%@*.
                  .=@@%@%##@@+----%@*.              .*@%----+@@##%@%@@=.
                .-@@%***@@%*#@@=----@#.           ..#@----=@@#*%@@**#@@@-.
              .:%@#@%*****@@@@*@%=--@#.            .#@--=%@*@@@@*****@@#@%:.
            .:#@%*#@%**###%@+...=%@%+.              .+%@%=. .+@%*****%@**#@#:.
          ..*@@@@@@@@@@@@@#.                                 ..#@@@@@@@@@@@@@*..
        ..+@@*****#@#*#@%.                                      .%@#**@#*****@@+.
    ..#%#@@*******#@%@%-..                                      ..-%@#@%*******@@#%#..
    +@#=%@#*******#@@-.         Once unlocked you will be         ..-@@@*******#@%=#@+
    %@=--=@@#****@@=..                able to run any               ..=@@****#@@=--=@%
    .-@@=--+@@#%@*..                  DRAKKAR command.               ..*@%#@@=--=@@-.
      .+@%---+@@:.                                                     .:@@+---%@+.
        .*@#--%@:                                                       .-@%--#@*..
          .#@@*.                                                         ..*@@#.

    """
    print(ascii_swords)

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

def is_snakemake_locked(workdir: str) -> bool:
    locks_dir = os.path.join(workdir, ".snakemake", "locks")
    return os.path.isdir(locks_dir) and len(os.listdir(locks_dir)) > 0

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

# Create dictionary of bin names and paths from the input file
def file_bins_to_json(paths_file=None, output=False):
    fasta_dict = {}

    if not os.path.isfile(paths_file):
            raise FileNotFoundError(f"Bin file not found: {paths_file}")

    # Read the paths file
    with open(paths_file, "r") as f:
        for line in f:
            full_path = line.strip()
            if not full_path:
                continue

            # Extract filename without path and extension
            filename = os.path.splitext(os.path.basename(full_path))[0]

            # Store in dictionary
            fasta_dict[filename] = full_path

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def path_bins_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"❌ Folder not found: {folder_path}")

    # Iterate over all files in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith((".fna", ".fa")):
            full_path = os.path.join(folder_path, file_name)
            file_id = os.path.splitext(file_name)[0]  # Remove extension
            fasta_dict[file_id] = full_path

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/bins_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

# Create dictionary of bin names and paths from the input file
def file_mags_to_json(paths_file=None, output=False):
    fasta_dict = {}

    # Read the paths file
    with open(paths_file, "r") as f:
        for line in f:
            full_path = line.strip()
            if not full_path:
                continue

            # Extract filename without path and extension
            filename = os.path.splitext(os.path.basename(full_path))[0]

            # Store in dictionary
            fasta_dict[filename] = full_path

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/mags_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

#updated to account for compressed genomes
def path_mags_to_json(folder_path=None, output=False):
    fasta_dict = {}

    # Ensure folder exists
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"❌ Folder not found: {folder_path}")

    # Compile a regex that matches .fa/.fna/.fasta, optionally followed by .gz
    FASTA_RE = re.compile(r'\.(?:fa|fna|fasta)(?:\.gz)?$', re.IGNORECASE)

    # Iterate over all files in the folder
    fasta_dict = {}
    for fname in os.listdir(folder_path):
        if FASTA_RE.search(fname):
            full_path = os.path.join(folder_path, fname)
            sample_id = FASTA_RE.sub('', fname)
            fasta_dict[sample_id] = full_path

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/mags_to_files.json", "w") as f:
        json.dump(fasta_dict, f, indent=4)

def microdiversity_selection_to_json(coverage, mincov=0.5, minsamp=10):
    df = pd.read_csv(csv_input, sep='\t', index_col=0)
    genome_sample_dict = {}
    for genome, row in df.iterrows():
        selected = row[row > mincov].index.tolist()
        if len(selected) >= minsamp:
            genome_sample_dict[genome] = selected

    with open(f"{output}/data/microdiversity_selection.json", "w") as f:
        json.dump(genome_sample_dict, f, indent=4)

def preprocessing_summary(summary_table, bar_width=50):
    """
    Prints a single horizontal stacked barplot representing the AVERAGE percentage of
    bases_discarded, bases_host, and bases_metagenomic across all samples.

    Uses:
    - '░' (light block) for discarded bases
    - '▒' (medium block) for host bases
    - '█' (full block) for metagenomic bases
    - '│' (vertical separator) between categories
    """

    df = pd.read_csv(summary_table, sep="\t")

    # Compute total sums
    total_discarded = df["bases_discarded"].sum()
    total_host = df["bases_host"].sum()
    total_metagenomic = df["bases_metagenomic"].sum()

    total_bases = total_discarded + total_host + total_metagenomic
    if total_bases == 0:
        print("No data available")
        return

    # Compute AVERAGE percentages
    pct_discarded = (total_discarded / total_bases) * 100
    pct_host = (total_host / total_bases) * 100
    pct_metagenomic = (total_metagenomic / total_bases) * 100

    # Compute character count for each section
    discarded_chars = round((pct_discarded / 100) * bar_width)
    host_chars = round((pct_host / 100) * bar_width)
    metagenomic_chars = bar_width - discarded_chars - host_chars  # Ensure total width matches

    # Construct the visual bar
    bar = "░" * discarded_chars + "│" + "▒" * host_chars + "│" + "█" * metagenomic_chars

    # Print the averages and the stacked barplot
    print(f"░ Discarded: {pct_discarded:.1f}%, ▒ Host: {pct_host:.1f}%, █ Metagenomic: {pct_metagenomic:.1f}%")
    print(f"│{bar}│")


def file_transcriptome_to_json(infofile, output):
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
    with open(f"{output}/data/transcriptome_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f)

    with open(f"{output}/data/transcriptome_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f)

def argument_transcriptome_to_json(argument, output):
    # Define the directory containing the raw reads
    TRANSCRIPTOME_DIR = Path(argument).resolve()

    # Initialize dictionaries
    TRANSCRIPTOME_TO_READS1 = defaultdict(list)
    TRANSCRIPTOME_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(TRANSCRIPTOME_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(TRANSCRIPTOME_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    TRANSCRIPTOME_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    TRANSCRIPTOME_TO_READS2[sample_name].append(full_path)

    # Convert defaultdict to standard dict (optional)
    TRANSCRIPTOME_TO_READS1 = dict(TRANSCRIPTOME_TO_READS1)
    TRANSCRIPTOME_TO_READS2 = dict(TRANSCRIPTOME_TO_READS2)

    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/transcriptome_to_reads1.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS1, f)

    with open(f"{output}/data/transcriptome_to_reads2.json", "w") as f:
        json.dump(TRANSCRIPTOME_TO_READS2, f)