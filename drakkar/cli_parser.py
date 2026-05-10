import argparse
import os

from drakkar import __version__
from drakkar.cli_context import DEFAULT_CATALOGING_BINNERS
from drakkar.cli_help import RichArgumentParser, _set_help_metadata
from drakkar.cli_validation import (
    add_benchmark_argument,
    add_resource_multiplier_arguments,
    positive_int,
)
from drakkar.database_registry import MANAGED_DATABASES


def build_parser():
    parser = RichArgumentParser(
        prog="drakkar",
        description="Drakkar: A Snakemake-based workflow for sequencing analysis",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--version", action="version", version=f"drakkar {__version__}")
    subparsers = parser.add_subparsers(dest="command", help="Available workflows")
    
    # Define subcommands for each workflow
    subparser_complete = subparsers.add_parser("complete", help="Run the end-to-end workflow from raw reads to catalog, profiling, and annotation outputs")
    subparser_complete.add_argument("-i", "--input", required=False, help="Input directory")
    subparser_complete.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_complete.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    complete_reference_group = subparser_complete.add_mutually_exclusive_group()
    complete_reference_group.add_argument("-r", "--reference", required=False, help="Reference host genome FASTA")
    complete_reference_group.add_argument("-x", "--reference-index", required=False, help="Tarball containing a reference FASTA and Bowtie2 index files")
    subparser_complete.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_complete.add_argument(
        "-b",
        "--binners",
        required=False,
        default=DEFAULT_CATALOGING_BINNERS,
        help=(
            "Comma-separated list of binners to use for binning. "
            "Options: metabat, maxbin, semibin, comebin. Default: all"
        ),
    )
    subparser_complete.add_argument("-t", "--type", required=False, default="genomes", help="Either genomes or pangenomes profiling type. Default: genomes")
    subparser_complete.add_argument(
        "--annotation-type",
        dest="annotation_type",
        required=False,
        default="taxonomy,function",
        help=(
            "Comma-separated annotation targets. Options: taxonomy, function, genes, "
            "kegg, cazy, pfam, virulence (vfdb), amr, signalp, dbcan, antismash, "
            "defense, mobile (genomad), network. Default: taxonomy,function"
        ),
    )
    subparser_complete.add_argument(
        "--gtdb-version",
        dest="gtdb_version",
        required=False,
        help="GTDB release number for taxonomy annotation, using GTDB_DB_<version> from config.yaml. Default: GTDB_DB.",
    )
    subparser_complete.add_argument("-c", "--multicoverage", action="store_true", help="Map samples sharing the same coverage group to each other's individual assemblies")
    subparser_complete.add_argument("--fraction", required=False, action='store_true', help="Calculate microbial fraction using singlem")
    subparser_complete.add_argument("--nonpareil", required=False, action='store_true', help="Estimate metagenomic coverage and diversity using Nonpareil during preprocessing")
    subparser_complete.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_complete.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_complete.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_complete.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_complete)
    add_resource_multiplier_arguments(subparser_complete)
    
    subparser_preprocessing = subparsers.add_parser("preprocessing", help="Quality-filter reads, optionally remove host sequences, and prepare cleaned datasets for downstream analysis")
    subparser_preprocessing.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_preprocessing.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_preprocessing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    preprocessing_reference_group = subparser_preprocessing.add_mutually_exclusive_group()
    preprocessing_reference_group.add_argument("-r", "--reference", required=False, help="Reference host genome FASTA")
    preprocessing_reference_group.add_argument("-x", "--reference-index", required=False, help="Tarball containing a reference FASTA and Bowtie2 index files")
    subparser_preprocessing.add_argument("--fraction", required=False, action='store_true', help="Calculate microbial fraction using singlem")
    subparser_preprocessing.add_argument("--nonpareil", required=False, action='store_true', help="Estimate metagenomic coverage and diversity using Nonpareil")
    subparser_preprocessing.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_preprocessing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_preprocessing.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_preprocessing)
    add_resource_multiplier_arguments(subparser_preprocessing)
    
    subparser_cataloging = subparsers.add_parser("cataloging", help="Assemble reads, bin genomes, and build the MAG catalog used by downstream workflows")
    subparser_cataloging.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_cataloging.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_cataloging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_cataloging.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_cataloging.add_argument(
        "-b",
        "--binners",
        required=False,
        default=DEFAULT_CATALOGING_BINNERS,
        help=(
            "Comma-separated list of binners to use for binning. "
            "Options: metabat, maxbin, semibin, comebin. Default: all"
        ),
    )
    subparser_cataloging.add_argument("-c", "--multicoverage", action="store_true", help="Map samples sharing the same coverage group to each other's individual assemblies")
    subparser_cataloging.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_cataloging.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_cataloging.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_cataloging)
    add_resource_multiplier_arguments(subparser_cataloging)
    
    subparser_profiling = subparsers.add_parser("profiling", help="Dereplicate MAGs and quantify genome or pangenome abundance across metagenomic samples")
    subparser_profiling.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_profiling.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_profiling.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_profiling.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_profiling.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_profiling.add_argument("-t", "--type", required=False, default="genomes", help="Either genomes or pangenomes profiling type. Default: genomes")
    subparser_profiling.add_argument("-f", "--fraction", required=False, action='store_true', help="Calculate microbial fraction using singlem")
    subparser_profiling.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_profiling.add_argument("-n", "--ignore_quality", action="store_true", help="Pass --ignoreGenomeQuality to dRep during profiling")
    subparser_profiling.add_argument("-q", "--quality", type=str, help="CSV/TSV with genome, completeness, contamination columns")
    subparser_profiling.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_profiling.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_profiling.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_profiling)
    add_resource_multiplier_arguments(subparser_profiling)
    
    subparser_dereplicating = subparsers.add_parser("dereplicating", help="Dereplicate genome collections and export representative MAG FASTA files without mapping reads")
    subparser_dereplicating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_dereplicating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_dereplicating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_dereplicating.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_dereplicating.add_argument("-n", "--ignore_quality", action="store_true", help="Pass --ignoreGenomeQuality to dRep during dereplication")
    subparser_dereplicating.add_argument("-q", "--quality", type=str, help="CSV/TSV with genome, completeness, contamination columns")
    subparser_dereplicating.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_dereplicating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_dereplicating.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_dereplicating)
    add_resource_multiplier_arguments(subparser_dereplicating)
    
    subparser_annotating = subparsers.add_parser("annotating", help="Assign taxonomy and functional annotations to MAGs, genes, and derived feature tables")
    subparser_annotating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa, .fna or .fasta, optionally including .gz) are stored")
    subparser_annotating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_annotating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_annotating.add_argument(
        "--annotation-type",
        dest="annotation_type",
        required=False,
        default="taxonomy,function",
        help=(
            "Comma-separated annotation targets. Options: taxonomy, function, genes, "
            "kegg, cazy, pfam, virulence (vfdb), amr, signalp, dbcan, antismash, "
            "defense, mobile (genomad), network. Default: taxonomy,function"
        ),
    )
    subparser_annotating.add_argument(
        "--gtdb-version",
        dest="gtdb_version",
        required=False,
        help="GTDB release number for taxonomy annotation, using GTDB_DB_<version> from config.yaml. Default: GTDB_DB.",
    )
    subparser_annotating.add_argument("-e", "--env_path", type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_annotating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_annotating.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_annotating)
    add_resource_multiplier_arguments(subparser_annotating)
    
    subparser_inspecting = subparsers.add_parser("inspecting", help="Combine bins and coverage or mapping inputs into inspection-ready summaries for follow-up exploration")
    subparser_inspecting.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_inspecting.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_inspecting.add_argument("-m", "--mapping_dir", required=False, help="Directory containing the mapping (.bam) files")
    subparser_inspecting.add_argument("-c", "--cov_file", required=False, help="Tab-separated file containing mapping coverage per genome per smaple")
    subparser_inspecting.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_inspecting.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_inspecting.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_inspecting.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_inspecting)
    add_resource_multiplier_arguments(subparser_inspecting)
    
    subparser_expressing = subparsers.add_parser("expressing", help="Quantify microbial gene expression from reads against MAG-derived gene predictions")
    subparser_expressing.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_expressing.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_expressing.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_expressing.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_expressing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_expressing.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_expressing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_expressing.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_benchmark_argument(subparser_expressing)
    add_resource_multiplier_arguments(subparser_expressing)
    
    database_parent = RichArgumentParser(add_help=False)
    database_parent.add_argument("--directory", required=True, help="Base directory where the database release directory will be created")
    database_parent.add_argument("--version", required=False, help="Release folder name to create inside --directory (optional for vfdb; defaults to the UTC download date)")
    database_parent.add_argument("--download-runtime", type=int, default=120, help="Runtime in minutes for the database download/preparation rule (default: 120)")
    database_parent.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    database_parent.add_argument("--set-default", action="store_true", help="Update config.yaml to use this installed database release by default")
    database_parent.add_argument("-e", "--env_path", type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    database_parent.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    add_benchmark_argument(database_parent)
    add_resource_multiplier_arguments(database_parent)
    
    subparser_database = subparsers.add_parser("database", help="Install or update managed annotation database releases used by Drakkar workflows")
    database_subparsers = subparser_database.add_subparsers(dest="database_name", help="Managed databases")
    database_subparsers.required = True
    
    database_kegg = database_subparsers.add_parser("kegg", parents=[database_parent], help="Install or update the KEGG/KOfam database", aliases=["kofams"])
    database_cazy = database_subparsers.add_parser("cazy", parents=[database_parent], help="Install or update the CAZy database")
    database_pfam = database_subparsers.add_parser("pfam", parents=[database_parent], help="Install or update the PFAM database")
    database_vfdb = database_subparsers.add_parser("vfdb", parents=[database_parent], help="Install or update the VFDB database")
    database_amr = database_subparsers.add_parser("amr", parents=[database_parent], help="Install or update the AMR database")
    
    subparser_environments = subparsers.add_parser("environments", help="Pre-build Drakkar conda environments before launching workflow jobs")
    subparser_environments.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_environments.add_argument("--profile", default="local", choices=["local", "slurm"])
    add_benchmark_argument(subparser_environments)
    add_resource_multiplier_arguments(subparser_environments)
    
    subparser_unlock = subparsers.add_parser("unlock", help="Remove a stale Snakemake lock from a Drakkar output directory")
    subparser_unlock.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_unlock.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_unlock.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    
    subparser_update = subparsers.add_parser("update", help="Reinstall Drakkar from GitHub in the current Python environment")
    subparser_update.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_update.add_argument("--skip-deps", action="store_true", help="Skip reinstalling Python package dependencies during the update")
    
    subparser_transfer = subparsers.add_parser("transfer", help="Transfer selected Drakkar outputs to a remote SFTP destination while preserving folder structure")
    subparser_transfer.add_argument("--host", help="SFTP host")
    subparser_transfer.add_argument("--user", help="SFTP user")
    subparser_transfer.add_argument("--port", type=int, default=22, help="SFTP port")
    subparser_transfer.add_argument("-i", "--identity", help="Path to SSH private key")
    subparser_transfer.add_argument("-l", "--local-dir", required=False, default=os.getcwd(), help="Local output directory. Default is the directory from which drakkar is called.")
    subparser_transfer.add_argument("-r", "--remote-dir", required=True, help="Remote base directory for transfer")
    subparser_transfer.add_argument("--erda", action="store_true", help="Use ERDA SFTP defaults")
    subparser_transfer.add_argument("--all", action="store_true", help="Transfer the entire output folder")
    subparser_transfer.add_argument("--data", action="store_true", help="Transfer the output folder excluding .snakemake")
    subparser_transfer.add_argument("--results", action="store_true", help="Transfer configured results files")
    subparser_transfer.add_argument("-a", "--annotations", action="store_true", help="Transfer annotation outputs")
    subparser_transfer.add_argument("-m", "--mags", action="store_true", help="Transfer dereplicated MAGs")
    subparser_transfer.add_argument("-p", "--profile", action="store_true", help="Transfer profiling outputs")
    subparser_transfer.add_argument("-e", "--expression", action="store_true", help="Transfer expression outputs")
    subparser_transfer.add_argument("-b", "--bins", action="store_true", help="Transfer cataloging bins")
    subparser_transfer.add_argument("-v", "--verbose", action="store_true", help="Log each transfer on screen")
    
    subparser_config = subparsers.add_parser("config", help="View or edit the installed workflow/config.yaml used by this Drakkar environment")
    config_actions = subparser_config.add_mutually_exclusive_group(required=True)
    config_actions.add_argument("--view", action="store_true", help="Print workflow/config.yaml")
    config_actions.add_argument("--edit", action="store_true", help="Open workflow/config.yaml in a terminal editor")
    
    subparser_logging = subparsers.add_parser("logging", help="Inspect run metadata, progress summaries, and Snakemake logs to troubleshoot workflows")
    subparser_logging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_logging.add_argument("--run", required=False, help="Specific run ID (YYYYMMDD-HHMMSS) or drakkar_<run_id>.yaml file name")
    subparser_logging.add_argument("--tail", required=False, type=positive_int, default=50, help="Number of log lines to show when --excerpt is used and no failure excerpt is found. Default: 50")
    subparser_logging.add_argument("--summary", action="store_true", help="Print the parsed workflow summary without showing log content")
    subparser_logging.add_argument("--excerpt", action="store_true", help="Print the most recent failure excerpt, or the last --tail lines if no excerpt is found")
    subparser_logging.add_argument("--full", action="store_true", help="Print the full Snakemake log")
    subparser_logging.add_argument("--paths", action="store_true", help="List relevant metadata and log paths")
    subparser_logging.add_argument("--list", action="store_true", help="List available workflow runs in the output directory")
    
    parser.description = "Genome-resolved metagenomics workflows, database setup, and run management."
    _set_help_metadata(
        parser,
        category="Workflow launcher and operations hub",
        examples=[
            "drakkar complete -f input_info.tsv -o drakkar_output",
            "drakkar preprocessing -i reads/ -o drakkar_output",
            "drakkar logging -o drakkar_output --summary",
        ],
        command_groups=[
            ("Data Generation and Analysis", ["complete", "preprocessing", "cataloging", "profiling", "dereplicating", "annotating", "inspecting", "expressing"]),
            ("Operations and Management", ["database", "environments", "logging", "transfer", "config", "unlock", "update"]),
        ],
        sections=[
            ("General", ["help", "--version"]),
        ],
    )
    
    subparser_complete.description = "Run preprocessing, cataloging, profiling, and annotation as a single end-to-end workflow."
    _set_help_metadata(
        subparser_complete,
        category="End-to-end workflow",
        examples=[
            "drakkar complete -f input_info.tsv -o drakkar_output",
            "drakkar complete -i reads/ -o drakkar_output -m individual,all",
            "drakkar complete -f input_info.tsv --annotation-type taxonomy,genes -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["input", "file", "reference", "reference_index"]),
            ("Workflow Scope", ["mode", "binners", "type", "annotation_type", "gtdb_version", "multicoverage", "fraction", "nonpareil", "ani"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_preprocessing.description = "Quality-filter reads, optionally remove host reads, and prepare datasets for downstream workflows."
    _set_help_metadata(
        subparser_preprocessing,
        category="Data generation workflow",
        examples=[
            "drakkar preprocessing -i reads/ -o drakkar_output",
            "drakkar preprocessing -f input_info.tsv -r host.fna -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["input", "file", "reference", "reference_index"]),
            ("Optional Analyses", ["fraction", "nonpareil"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_cataloging.description = "Assemble reads, bin genomes, and build the MAG catalog used by later workflows."
    _set_help_metadata(
        subparser_cataloging,
        category="Data generation workflow",
        examples=[
            "drakkar cataloging -i reads/ -o drakkar_output",
            "drakkar cataloging -f input_info.tsv -m individual,all -c -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["input", "file"]),
            ("Assembly Strategy", ["mode", "binners", "multicoverage"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_profiling.description = "Dereplicate MAGs and quantify genomes or pangenomes across metagenomic samples."
    _set_help_metadata(
        subparser_profiling,
        category="Analysis workflow",
        examples=[
            "drakkar profiling -b genomes/ -R input_info.tsv -o drakkar_output",
            "drakkar profiling -o drakkar_output -a 0.95",
            "drakkar profiling -B bins.txt -R input_info.tsv -q mag_qualities.tsv -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["bins_dir", "bins_file", "reads_dir", "reads_file"]),
            ("Analysis Settings", ["type", "fraction", "ani", "ignore_quality", "quality"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_dereplicating.description = "Run only MAG dereplication and export dereplicated genomes without read mapping."
    _set_help_metadata(
        subparser_dereplicating,
        category="Analysis workflow",
        examples=[
            "drakkar dereplicating -b genomes/ -o drakkar_output",
            "drakkar dereplicating -B bins.txt -q mag_qualities.csv -o drakkar_output",
        ],
        sections=[
            ("Input Genomes", ["bins_dir", "bins_file"]),
            ("Dereplication Settings", ["ani", "ignore_quality", "quality"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_annotating.description = "Annotate MAGs with taxonomy and selected functional modules."
    _set_help_metadata(
        subparser_annotating,
        category="Analysis workflow",
        examples=[
            "drakkar annotating -b genomes/ -o drakkar_output",
            "drakkar annotating -B bins.txt --annotation-type taxonomy,kegg,dbcan -o drakkar_output",
        ],
        sections=[
            ("Input Genomes", ["bins_dir", "bins_file"]),
            ("Annotation Scope", ["annotation_type", "gtdb_version"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_inspecting.description = "Combine bins and mapping coverage into inspection-ready summaries."
    _set_help_metadata(
        subparser_inspecting,
        category="Analysis workflow",
        examples=[
            "drakkar inspecting -b genomes/ -m mapping/ -o drakkar_output",
            "drakkar inspecting -B bins.txt -c coverage.tsv -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["bins_dir", "bins_file", "mapping_dir", "cov_file"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_expressing.description = "Quantify microbial gene expression from reads and MAG references."
    _set_help_metadata(
        subparser_expressing,
        category="Analysis workflow",
        examples=[
            "drakkar expressing -b genomes/ -R input_info.tsv -o drakkar_output",
            "drakkar expressing -B bins.txt -r reads/ -o drakkar_output",
        ],
        sections=[
            ("Input Sources", ["bins_dir", "bins_file", "reads_dir", "reads_file"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_database.description = "Install or update one managed annotation database release and optionally make it the default in config.yaml."
    _set_help_metadata(
        subparser_database,
        category="Operations and management",
        examples=[
            "drakkar database kegg --directory /db/kofams --version 2026-02-01",
            "drakkar database vfdb --directory /db/vfdb --set-default",
        ],
        command_groups=[
            ("Managed Databases", ["kegg", "cazy", "pfam", "vfdb", "amr"]),
        ],
        sections=[
            ("General", ["help"]),
        ],
    )
    
    database_descriptions = {
        database_kegg: (
            "Install a versioned KEGG/KOfam profile release and prepare the pressed HMM database.",
            [
                "drakkar database kegg --directory /db/kofams --version 2026-02-01",
                "drakkar database kegg --directory /db/kofams --version 2026-02-01 --set-default",
            ],
        ),
        database_cazy: (
            "Install a dbCAN HMM release for CAZy-style annotation and prepare the pressed HMM database.",
            [
                "drakkar database cazy --directory /db/cazy --version V14",
                "drakkar database cazy --directory /db/cazy --version V14 --set-default",
            ],
        ),
        database_pfam: (
            "Install a Pfam release and prepare the pressed HMM database used by annotation rules.",
            [
                "drakkar database pfam --directory /db/pfam --version Pfam37.4",
            ],
        ),
        database_vfdb: (
            "Install the latest VFDB protein set and store it under a download-date release folder.",
            [
                "drakkar database vfdb --directory /db/vfdb",
                "drakkar database vfdb --directory /db/vfdb --set-default",
            ],
        ),
        database_amr: (
            "Install a versioned NCBIfam-AMRFinder release used for antimicrobial resistance annotation.",
            [
                "drakkar database amr --directory /db/amr --version 2025-07-16.1",
            ],
        ),
    }
    for db_parser, (description, examples) in database_descriptions.items():
        db_parser.description = description
        _set_help_metadata(
            db_parser,
            category="Database management",
            examples=examples,
            sections=[
                ("Release Settings", ["directory", "version", "download_runtime", "set_default"]),
                ("Run Configuration", ["env_path", "profile", "overwrite", "skip_benchmark"]),
                ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
            ],
        )
    
    subparser_environments.description = "Pre-build workflow conda environments before launching production runs."
    _set_help_metadata(
        subparser_environments,
        category="Operations and management",
        examples=[
            "drakkar environments --profile local",
            "drakkar environments -e /shared/drakkar_envs --profile slurm",
        ],
        sections=[
            ("Environment Setup", ["env_path", "profile", "skip_benchmark"]),
            ("Resource Scaling", ["memory_multiplier", "time_multiplier"]),
        ],
    )
    
    subparser_unlock.description = "Clear a Snakemake lock from an output directory when you are sure no workflow is running."
    _set_help_metadata(
        subparser_unlock,
        category="Operations and management",
        examples=[
            "drakkar unlock -o drakkar_output",
        ],
        sections=[
            ("Target Directory", ["output"]),
            ("Execution", ["env_path", "profile"]),
        ],
    )
    
    subparser_update.description = "Reinstall the current Drakkar package directly from GitHub in the active Python environment."
    _set_help_metadata(
        subparser_update,
        category="Operations and management",
        examples=[
            "drakkar update",
            "drakkar update --skip-deps",
        ],
        sections=[
            ("Execution", ["env_path", "skip_deps"]),
        ],
    )
    
    subparser_transfer.description = "Copy selected Drakkar outputs to an SFTP destination while preserving the project folder structure."
    _set_help_metadata(
        subparser_transfer,
        category="Operations and management",
        examples=[
            "drakkar transfer -l drakkar_output -r remote/project --results --annotations",
            "drakkar transfer --erda -l drakkar_output -r holocamp/project --data -v",
        ],
        sections=[
            ("Connection", ["host", "user", "port", "identity", "erda"]),
            ("Locations", ["local_dir", "remote_dir"]),
            ("What to Transfer", ["all", "data", "results", "annotations", "mags", "profile", "expression", "bins"]),
            ("Display", ["verbose"]),
        ],
    )
    
    subparser_config.description = "View or edit the installed workflow/config.yaml used by the current Drakkar installation."
    _set_help_metadata(
        subparser_config,
        category="Operations and management",
        examples=[
            "drakkar config --view",
            "drakkar config --edit",
        ],
        sections=[
            ("Actions", ["view", "edit"]),
        ],
    )
    
    subparser_logging.description = "Inspect run metadata and Snakemake logs to troubleshoot failed, interrupted, or incomplete workflows."
    _set_help_metadata(
        subparser_logging,
        category="Operations and management",
        examples=[
            "drakkar logging -o drakkar_output",
            "drakkar logging -o drakkar_output --summary",
            "drakkar logging -o drakkar_output --excerpt",
            "drakkar logging -o drakkar_output --run 20260503-101530 --paths",
        ],
        sections=[
            ("Target Run", ["output", "run"]),
            ("Display Options", ["summary", "excerpt", "tail", "full", "paths", "list"]),
        ],
    )
    
    return parser
