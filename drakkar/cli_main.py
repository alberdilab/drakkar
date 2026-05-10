import json
import os
from pathlib import Path

from drakkar.cli_context import (
    ERROR,
    INFO,
    READ_ONLY_COMMANDS,
    RESET,
    WORKFLOW_RUN_COMMANDS,
    config_vars,
)
from drakkar.cli_parser import build_parser
from drakkar.cli_validation import (
    normalize_annotation_type,
    normalize_cataloging_binners,
    validate_download_runtime,
    validate_gtdb_version,
    validate_managed_database_version,
)
from drakkar.config_commands import edit_config, set_default_database_path, view_config
from drakkar.database_registry import MANAGED_DATABASES, database_release_dir, normalize_managed_database_name
from drakkar.output import print, section
from drakkar.output_paths import prepare_output_directory, validate_path
from drakkar.quality import validate_and_write_quality_file
from drakkar.run_logs import run_logging
from drakkar.run_metadata import write_launch_metadata
from drakkar.system_checks import check_screen_session
from drakkar.transfer import run_sftp_transfer
from drakkar.update_command import run_update
from drakkar.utils import (
    argument_preprocessed_to_json,
    argument_references_to_json,
    argument_samples_to_json,
    argument_transcriptome_to_json,
    check_reference_columns,
    display_drakkar,
    file_assemblies_to_json,
    file_bins_to_json,
    file_coverages_to_json,
    file_mags_to_json,
    file_preprocessed_to_json,
    file_references_to_json,
    file_samples_to_json,
    file_transcriptome_to_json,
    path_bins_to_json,
    path_mags_to_json,
)
from drakkar.workflow_launcher import (
    run_snakemake_annotating,
    run_snakemake_cataloging,
    run_snakemake_database,
    run_snakemake_dereplicating,
    run_snakemake_environments,
    run_snakemake_expressing,
    run_snakemake_inspecting,
    run_snakemake_preprocessing,
    run_snakemake_profiling,
    run_unlock,
)

def main():
    parser = build_parser()
    args = parser.parse_args()

    # Display ASCII logo before running any command or showing help
    display_drakkar()

    # Check screen session
    if args.command not in READ_ONLY_COMMANDS:
        check_screen_session()

    path_checks = [
        (getattr(args, "input", None), "Input", True),
        (getattr(args, "file", None), "Sample detail file", False),
        (getattr(args, "bins_dir", None), "Bins directory", True),
        (getattr(args, "bins_file", None), "Bins file", False, True),
        (getattr(args, "reads_dir", None), "Reads directory", True),
        (getattr(args, "reads_file", None), "Reads file", False),
        (getattr(args, "reference", None), "Reference", False, True),
        (getattr(args, "reference_index", None), "Reference index tarball", False, True),
        (getattr(args, "cov_file", None), "Coverage file", False),
    ]
    for path_check in path_checks:
        path_value, label, expect_dir, *options = path_check
        allow_url = options[0] if options else False
        if not validate_path(path_value, label, expect_dir, allow_url=allow_url):
            return

    ###
    # Declare environment directory
    ###

    env_path = getattr(args, "env_path", None)
    if env_path:
        env_path = env_path
    else:
        if args.command != "update":
            env_path = config_vars['ENVIRONMENTS_DIR']

    if args.command in ("annotating", "complete"):
        normalized_annotation_type = normalize_annotation_type(args.annotation_type)
        if not normalized_annotation_type:
            return
        args.annotation_type = normalized_annotation_type
        normalized_gtdb_version = validate_gtdb_version(getattr(args, "gtdb_version", None))
        if getattr(args, "gtdb_version", None) and not normalized_gtdb_version:
            return
        args.gtdb_version = normalized_gtdb_version

    if args.command in ("cataloging", "complete"):
        normalized_binners = normalize_cataloging_binners(getattr(args, "binners", None))
        if not normalized_binners:
            return
        args.binners = normalized_binners

    if args.command == "database":
        normalized_database_name = normalize_managed_database_name(getattr(args, "database_name", None))
        if not normalized_database_name:
            print(f"{ERROR}ERROR:{RESET} Supported database commands are: {', '.join(MANAGED_DATABASES)}")
            return
        version_was_provided = bool(getattr(args, "version", None))
        normalized_database_version = validate_managed_database_version(normalized_database_name, getattr(args, "version", None))
        if not normalized_database_version:
            return
        normalized_download_runtime = validate_download_runtime(getattr(args, "download_runtime", None))
        if normalized_download_runtime is None:
            return
        if Path(args.directory).exists() and Path(args.directory).is_file():
            print(f"{ERROR}ERROR:{RESET} --directory must be a directory path, not a file: {args.directory}")
            return
        args.database_name = normalized_database_name
        args.version = normalized_database_version
        args.download_runtime = normalized_download_runtime
        if args.database_name == "vfdb" and not version_was_provided:
            print(f"{INFO}INFO:{RESET} No --version provided for vfdb; using download date {args.version} (UTC).")

    overwrite_capable_commands = {
        "complete", "preprocessing", "cataloging", "profiling", "dereplicating",
        "annotating", "inspecting", "expressing", "database",
    }

    if args.command == "transfer":
        output_dir = getattr(args, "local_dir", os.getcwd())
    elif args.command == "database":
        output_dir = database_release_dir(args.database_name, args.directory, args.version)
    else:
        output_dir = getattr(args, "output", os.getcwd())

    if args.command in overwrite_capable_commands:
        if not prepare_output_directory(output_dir, overwrite=getattr(args, "overwrite", False)):
            return
    run_info = None
    if args.command in WORKFLOW_RUN_COMMANDS:
        run_info = write_launch_metadata(args, output_dir, env_path=locals().get("env_path"))
        if not run_info:
            return

    ###
    # Unlock, update or create environments
    ###

    if args.command == "unlock":
        section("UNLOCKING DRAKKAR DIRECTORY")
        run_unlock(args.command, args.output, args.profile)

    elif args.command == "config":
        if args.view:
            return view_config()
        if args.edit:
            return edit_config()
        return 0

    elif args.command == "logging":
        return run_logging(
            args.output,
            run_id=args.run,
            tail=args.tail,
            summary=args.summary,
            excerpt=args.excerpt,
            full=args.full,
            paths=args.paths,
            list_runs=args.list,
        )

    elif args.command == "transfer":
        section("TRANSFERRING DRAKKAR OUTPUTS")
        run_sftp_transfer(args)
        return

    elif args.command == "update":
        section("UPDATING DRAKKAR")
        return run_update(skip_deps=args.skip_deps)

    elif args.command == "environments":
        section("CREATING CONDA ENVIRONMENTS")
        run_snakemake_environments(
            args.command,
            env_path,
            args.profile,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )

    elif args.command == "database":
        section("UPDATING DRAKKAR DATABASE")
        release_dir = database_release_dir(args.database_name, args.directory, args.version).resolve()
        project_name = os.path.basename(os.path.normpath(release_dir))
        run_snakemake_database(
            "database",
            project_name,
            release_dir,
            env_path,
            args.profile,
            args.database_name,
            Path(args.directory).resolve(),
            args.version,
            args.download_runtime,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )
        if args.set_default:
            default_path = set_default_database_path(args.database_name, Path(args.directory).resolve(), args.version)
            print(f"{INFO}INFO:{RESET} Updated {MANAGED_DATABASES[args.database_name]['config_key']} in config.yaml to {default_path}")
        return

    else:            
        project_name = os.path.basename(os.path.normpath(args.output))

    ###
    # Preprocessing
    ###

    if args.command in ("preprocessing", "complete"):
        section("STARTING PREPROCESSING PIPELINE")

        # Generate raw data dictionaries

        if args.file and args.input:
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.file and not args.input:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_samples_to_json(args.file,args.output)

        elif args.input and not args.file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_samples_to_json(args.input,args.output)

        else:
            print(f"Please provide either an input directory (-i) or a sample info file (-f)")
            return

        # Generate reference genome dictionaries
        reference_argument = getattr(args, "reference", None) or getattr(args, "reference_index", None)
        reference_source_label = "reference index tarball" if getattr(args, "reference_index", None) else "reference genome file"

        if args.file and reference_argument:
            if check_reference_columns(args.file):
                print(f"")
                print(f"Both sample info file and {reference_source_label} were provided.")
                print(f"DRAKKAR will continue with the information provided in the sample info file.")
                file_references_to_json(args.file,args.output)
                REFERENCE = True
            else:
                argument_references_to_json(reference_argument,f"{args.output}/data/sample_to_reads1.json",args.output)
                REFERENCE = True

        elif args.file and not reference_argument:
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

        elif reference_argument and not args.file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the {reference_source_label}.")
            argument_references_to_json(reference_argument,f"{args.output}/data/sample_to_reads1.json",args.output)
            REFERENCE = True

        else:
            print(f"")
            print(f"Running DRAKKAR without mapping against a reference genome")
            REFERENCE = False

        run_snakemake_preprocessing(
            "preprocessing",
            project_name,
            Path(args.output).resolve(),
            REFERENCE,
            env_path,
            args.profile,
            args.fraction,
            args.nonpareil,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )

    ###
    # Cataloging
    ###

    if args.command in ("cataloging", "complete"):
        section("STARTING CATALOGING PIPELINE")

        # Generate cataloging data dictionaries

        if args.file and args.input:
            print(f"Both sample info file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        elif args.file and not args.input:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_preprocessed_to_json(args.file,args.output)
        elif args.input and not args.file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            argument_preprocessed_to_json(args.input,args.output)
        else:
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

        with open(f"{args.output}/data/assembly_to_samples.json", "r") as f:
            _assembly_check = json.load(f)
        if not _assembly_check and not args.mode:
            print(f"")
            print(f"No assembly information found in the sample info file and no assembly mode (-m) was provided.")
            print(f"    In consequence, DRAKKAR will run individual assemblies.")
            INDIVIDUAL_MODE=True
            file_assemblies_to_json(args.file,samples,INDIVIDUAL_MODE,ALL_MODE,args.output)

        if args.multicoverage:
            with open(f"{args.output}/data/assembly_to_samples.json", "r") as f:
                assembly_to_samples = json.load(f)
            if any(len(group_samples) > 1 for group_samples in assembly_to_samples.values()):
                print(f"{ERROR}ERROR:{RESET} --multicoverage is not compatible with co-assemblies.")
                print("Co-assemblies already use coverage from the samples used to build the assembly.")
                return
            has_coverage = file_coverages_to_json(args.file, samples, args.output)
            if has_coverage:
                print("Coverage groups loaded from the sample info file.")
            else:
                print("No coverage information provided; all samples will be mapped against all assemblies.")

        run_snakemake_cataloging(
            "cataloging",
            project_name,
            Path(args.output).resolve(),
            env_path,
            args.profile,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
            binners=args.binners,
        )

    ###
    # Profiling
    ###

    if args.command in ("profiling", "complete"):
        section("STARTING PROFILING PIPELINE")

        # Prepare bin dictionaries

        if args.command in ("profiling"):
            bins_dir = args.bins_dir
            bins_file = args.bins_file

        if args.command in ("complete"):
            bins_dir = None
            bins_file = None

        if bins_dir and bins_file:
            print(f"Both bin path file and bin directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the bin info file.")
            print(bins_file)
            file_bins_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"DRAKKAR will run with the bins in the bin directory.")
            path_bins_to_json(bins_dir,args.output)
        else:
            print(f"No bin information was provided. DRAKKAR will try to guess the location of the bins.")
            if os.path.exists(f"{args.output}/cataloging/final/all_bin_paths.txt"):
                file_bins_to_json(f"{args.output}/cataloging/final/all_bin_paths.txt",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an bin file (-B) or directory (-b).")
                return

        # Prepare read dictionaries

        if args.command in ("profiling"):
            reads_dir = args.reads_dir
            reads_file = args.reads_file

        if args.command in ("complete"):
            reads_dir = None
            reads_file = None

        if reads_dir and reads_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_preprocessed_to_json(reads_file,args.output)
        elif reads_file and not reads_dir:
            print(f"")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_preprocessed_to_json(reads_file,args.output)
        elif reads_dir and not reads_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reads from the input directory.")
            argument_preprocessed_to_json(reads_dir,args.output)
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

        quality_file = getattr(args, "quality", None)
        ignore_quality = args.ignore_quality
        if quality_file:
            if not validate_and_write_quality_file(quality_file, args.output):
                return
            ignore_quality = True
        run_snakemake_profiling(
            "profiling",
            project_name,
            args.type,
            args.output,
            env_path,
            args.profile,
            args.fraction,
            args.ani,
            ignore_quality,
            quality_file,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )

    ###
    # Dereplicating
    ###

    if args.command == "dereplicating":
        section("STARTING DEREPLICATING PIPELINE")

        bins_dir = args.bins_dir
        bins_file = args.bins_file

        if bins_dir and bins_file:
            print(f"Both bin path file and bin directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the bin info file.")
            file_bins_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"DRAKKAR will run with the bins in the bin directory.")
            path_bins_to_json(bins_dir,args.output)
        else:
            print(f"No bin information was provided. DRAKKAR will try to guess the location of the bins.")
            if os.path.exists(f"{args.output}/cataloging/final/all_bin_paths.txt"):
                file_bins_to_json(f"{args.output}/cataloging/final/all_bin_paths.txt",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the cataloging module was run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an bin file (-B) or directory (-b).")
                return

        quality_file = getattr(args, "quality", None)
        ignore_quality = args.ignore_quality
        if quality_file:
            if not validate_and_write_quality_file(quality_file, args.output):
                return
            ignore_quality = True
        run_snakemake_dereplicating(
            "dereplicating",
            project_name,
            args.output,
            env_path,
            args.profile,
            args.ani,
            ignore_quality,
            quality_file,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )

    ###
    # Annotating
    ###

    if args.command in ("annotating", "complete"):
        section("STARTING ANNOTATING PIPELINE")

        if args.command in ("annotating"):
            bins_dir = args.bins_dir
            bins_file = args.bins_file

        if args.command in ("complete"):
            bins_dir = None
            bins_file = None

        # Prepare bin dictionaries
        if bins_dir and bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(bins_file,args.output)
        elif bins_file and not bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(bins_file,args.output)
        elif bins_dir and not bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        run_snakemake_annotating(
            "annotating",
            project_name,
            args.annotation_type,
            args.output,
            env_path,
            args.profile,
            getattr(args, "gtdb_version", None),
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )

    ###
    # Inspecting
    ###

    if args.command == "inspecting":
        section("STARTING INSPECTING PIPELINE")

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(args.bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return
            
        run_snakemake_inspecting(
            "inspecting",
            project_name,
            args.output,
            env_path,
            args.profile,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )


    ###
    # Expressing
    ###

    if args.command == "expressing":
        section("STARTING EXPRESSING PIPELINE")

        # Prepare bin dictionaries
        if args.bins_dir and args.bins_file:
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will continue with the information provided in the path file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_file and not args.bins_dir:
            print(f"DRAKKAR will run with the information provided in the sample info file.")
            file_mags_to_json(args.bins_file,args.output)
        elif args.bins_dir and not args.bins_file:
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will run with the files in the input directory.")
            path_mags_to_json(args.bins_dir,args.output)
        else:
            print(f"No input information was provided. DRAKKAR will try to guess the location of the MAGs.")
            if os.path.exists(f"{args.output}/profiling_genomes/drep/dereplicated_genomes"):
                path_mags_to_json(f"{args.output}/profiling_genomes/drep/dereplicated_genomes",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        # Prepare read dictionaries

        reads_dir = args.reads_dir
        reads_file = args.reads_file

        if reads_dir and reads_file:
            print(f"")
            print(f"Both bin path file and input directory were provided.")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_transcriptome_to_json(reads_file,args.output)
        elif reads_file and not reads_dir:
            print(f"")
            print(f"DRAKKAR will use the reads provided in the file.")
            file_transcriptome_to_json(reads_file,args.output)
        elif reads_dir and not reads_file:
            print(f"")
            print(f"No sample info file was provided.")
            print(f"DRAKKAR will use the reads from the input directory.")
            argument_transcriptome_to_json(reads_dir,args.output)
        else:
            print(f"")
            print(f"No input information was provided. DRAKKAR will try to guess the location of the reads.")
            if any(os.scandir(f"{args.output}/preprocessing/final")):
                argument_transcriptome_to_json(f"{args.output}/preprocessing/final",args.output)
            else:
                print(f"ERROR: No bin data was found in the output directory.")
                print(f"Make sure that the preprocessing and cataloging modules were run in this directory.")
                print(f"If you want to start from your own bin files, make sure to indicate an input file (-f) or directory (-i).")
                return

        run_snakemake_expressing(
            "expressing",
            project_name,
            args.output,
            env_path,
            args.profile,
            args.memory_multiplier,
            args.time_multiplier,
            run_info,
        )
