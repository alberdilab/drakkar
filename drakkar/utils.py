import os

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

def file_samples_to_json(infofile, output):
    # Load sample info file
    df = pd.read_csv(infofile, sep="\t")

    # Initialize dictionaries with lists
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Populate the dictionaries
    for _, row in df.iterrows():
        SAMPLE_TO_READS1[row["sample"]].append(row["rawreads1"])
        SAMPLE_TO_READS2[row["sample"]].append(row["rawreads2"])

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    # Print for verification
    print("SAMPLE_TO_READS1:", SAMPLE_TO_READS1)
    print("SAMPLE_TO_READS2:", SAMPLE_TO_READS2)

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
    with open(f"{args.output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    SAMPLE_TO_REFERENCE = dict(zip(df["sample"], df["reference_name"]))
    with open(f"{args.output}/data/sample_to_reference.json", "w") as f:
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
    REFERENCE_TO_FILE = {"reference": [f"2{argument}"]}
    os.makedirs(f"{output}/data", exist_ok=True)
    with open(f"{output}/data/reference_to_file.json", "w") as f:
        json.dump(REFERENCE_TO_FILE, f, indent=4)

    with open(f"{sample_to_reads}", "r") as f:
        SAMPLE_TO_READS = json.load(f)

    SAMPLE_TO_REFERENCE = {sample: "reference" for sample in SAMPLE_TO_READS.keys()}
    with open(f"{output}/data/sample_to_reference.json", "w") as f:
        json.dump(SAMPLE_TO_REFERENCE, f, indent=4)

#def argument_assemblies_to_json(argument):
