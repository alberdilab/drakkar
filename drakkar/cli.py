import argparse
import csv
import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import pandas as pd
import yaml

try:
    from importlib.metadata import PackageNotFoundError, version as get_distribution_version
except ImportError:  # pragma: no cover - Python < 3.8 fallback
    try:
        from importlib_metadata import PackageNotFoundError, version as get_distribution_version
    except ImportError:  # pragma: no cover - fallback if backport is absent
        PackageNotFoundError = Exception
        get_distribution_version = None

from drakkar import __version__
from drakkar import benchmark as _benchmark
from drakkar import cli_help as _cli_help
from drakkar import cli_main as _cli_main
from drakkar import cli_parser as _cli_parser
from drakkar import cli_validation as _cli_validation
from drakkar import config_commands as _config_commands
from drakkar import output_paths as _output_paths
from drakkar import run_logs as _run_logs
from drakkar import run_metadata as _run_metadata
from drakkar import update_command as _update_command
from drakkar import workflow_launcher as _workflow_launcher
from drakkar.cli_context import (
    CATALOGING_BINNER_ALIASES,
    CATALOGING_BINNER_ORDER,
    CONFIG_PATH,
    DEFAULT_CATALOGING_BINNERS,
    ERROR,
    HEADER1,
    INFO,
    PACKAGE_DIR,
    READ_ONLY_COMMANDS,
    RESET,
    WORKFLOW_RUN_COMMANDS,
    config_vars,
    load_config,
)
from drakkar.database_registry import (
    MANAGED_DATABASES,
    database_release_dir,
    normalize_managed_database_name,
)
from drakkar.output import get_console, print, prompt, section
from drakkar.utils import *

# Public Rich dependency aliases retained for tests and downstream monkeypatching.
rich_box = _cli_help.rich_box
RichGroup = _cli_help.RichGroup
RichPanel = _cli_help.RichPanel
RichTable = _cli_help.RichTable
RichText = _cli_help.RichText

# Direct re-exports that do not require compatibility shims.
_set_help_metadata = _cli_help._set_help_metadata
_clean_usage = _cli_help._clean_usage
_action_value = _cli_help._action_value
_expanded_action_help = _cli_help._expanded_action_help
_display_group_title = _cli_help._display_group_title
_find_action_by_key = _cli_help._find_action_by_key
_subparser_action = _cli_help._subparser_action

BENCHMARK_JOB_FIELDS = _benchmark.BENCHMARK_JOB_FIELDS
BENCHMARK_RULE_FIELDS = _benchmark.BENCHMARK_RULE_FIELDS
parse_int_or_none = _benchmark.parse_int_or_none
parse_float_or_none = _benchmark.parse_float_or_none
parse_resource_assignments = _benchmark.parse_resource_assignments
parse_slurm_memory_to_mb = _benchmark.parse_slurm_memory_to_mb
safe_ratio = _benchmark.safe_ratio
median_or_none = _benchmark.median_or_none
format_hours = _benchmark.format_hours
format_megabytes = _benchmark.format_megabytes
format_percent = _benchmark.format_percent
build_logical_job_key = _benchmark.build_logical_job_key
parse_snakemake_submitted_launches = _benchmark.parse_snakemake_submitted_launches
query_sacct_for_jobs = _benchmark.query_sacct_for_jobs
benchmark_job_row = _benchmark.benchmark_job_row
write_tsv = _benchmark.write_tsv
write_benchmark_summary_file = _benchmark.write_benchmark_summary_file
summarize_benchmark_rows = _benchmark.summarize_benchmark_rows
summarize_benchmark_rules = _benchmark.summarize_benchmark_rules
write_benchmark_tables = _benchmark.write_benchmark_tables
generate_benchmark_reports = _benchmark.generate_benchmark_reports
get_run_profile = _benchmark.get_run_profile
should_benchmark_run = _benchmark.should_benchmark_run
should_skip_benchmark = _benchmark.should_skip_benchmark

get_modules_to_run = _run_metadata.get_modules_to_run
build_snakemake_log_path = _run_metadata.build_snakemake_log_path
build_benchmark_paths = _run_metadata.build_benchmark_paths
load_metadata_file = _run_metadata.load_metadata_file
update_launch_metadata = _run_metadata.update_launch_metadata
finalize_launch_metadata = _run_metadata.finalize_launch_metadata
run_subprocess_with_logging = _run_metadata.run_subprocess_with_logging

workflow_run_sort_key = _run_logs.workflow_run_sort_key
discover_run_metadata = _run_logs.discover_run_metadata
resolve_run_metadata = _run_logs.resolve_run_metadata
discover_snakemake_fallback_logs = _run_logs.discover_snakemake_fallback_logs
tail_file = _run_logs.tail_file
extract_failure_excerpt = _run_logs.extract_failure_excerpt
classify_error_line = _run_logs.classify_error_line
summarize_snakemake_log = _run_logs.summarize_snakemake_log
print_snakemake_summary = _run_logs.print_snakemake_summary
print_benchmark_summary = _run_logs.print_benchmark_summary
print_logging_usage_guide = _run_logs.print_logging_usage_guide

from drakkar.quality import load_bins_map, normalize_genome_name, validate_and_write_quality_file
from drakkar.transfer import build_sftp_batch_commands, collect_transfer_files, list_files_recursive, run_sftp_transfer

_CONFIG_RESOLVE_EDITOR_COMMAND = _config_commands.resolve_editor_command
_CONFIG_REPLACE_CONFIG_VALUE = _config_commands.replace_config_value
_CONFIG_SET_DEFAULT_DATABASE_PATH = _config_commands.set_default_database_path
_CONFIG_VIEW_CONFIG = _config_commands.view_config
_CONFIG_EDIT_CONFIG = _config_commands.edit_config
_VALIDATION_DEFAULT_DATABASE_VERSION = _cli_validation.default_database_version
_UPDATE_GET_INSTALLED_VERSION = _update_command.get_installed_drakkar_version
_OUTPUT_OVERWRITE_OUTPUT_DIRECTORY = _output_paths.overwrite_output_directory
_OUTPUT_PROMPT_OVERWRITE_LOCKED_DIRECTORY = _output_paths.prompt_overwrite_locked_directory


def _is_facade_function(value, name):
    return (
        getattr(value, "__module__", None) == __name__
        and getattr(value, "__name__", None) == name
    )


def _sync_validation_dependencies():
    _cli_validation.config_vars = config_vars
    _cli_validation.print = print
    default_version_func = globals()["default_database_version"]
    if _is_facade_function(default_version_func, "default_database_version"):
        default_version_func = _VALIDATION_DEFAULT_DATABASE_VERSION
    _cli_validation.default_database_version = default_version_func


def _sync_config_dependencies():
    _config_commands.CONFIG_PATH = CONFIG_PATH
    _config_commands.print = print
    editor_func = globals()["resolve_editor_command"]
    if _is_facade_function(editor_func, "resolve_editor_command"):
        editor_func = _CONFIG_RESOLVE_EDITOR_COMMAND
    _config_commands.resolve_editor_command = editor_func


def _sync_help_dependencies():
    _cli_help.rich_box = rich_box
    _cli_help.RichGroup = RichGroup
    _cli_help.RichPanel = RichPanel
    _cli_help.RichTable = RichTable
    _cli_help.RichText = RichText
    _cli_help.print = print


def _sync_output_path_dependencies():
    _output_paths.print = print
    _output_paths.prompt = prompt
    _output_paths.is_snakemake_locked = globals()["is_snakemake_locked"]
    overwrite_func = globals()["overwrite_output_directory"]
    if _is_facade_function(overwrite_func, "overwrite_output_directory"):
        overwrite_func = _OUTPUT_OVERWRITE_OUTPUT_DIRECTORY
    prompt_func = globals()["prompt_overwrite_locked_directory"]
    if _is_facade_function(prompt_func, "prompt_overwrite_locked_directory"):
        prompt_func = _OUTPUT_PROMPT_OVERWRITE_LOCKED_DIRECTORY
    _output_paths.overwrite_output_directory = overwrite_func
    _output_paths.prompt_overwrite_locked_directory = prompt_func


def _sync_workflow_dependencies():
    _workflow_launcher.config_vars = config_vars
    _workflow_launcher.print = print
    _workflow_launcher.resource_config = globals()["resource_config"]
    _workflow_launcher.default_resource_args = globals()["default_resource_args"]
    subprocess_runner = globals()["run_subprocess_with_logging"]
    _workflow_launcher.run_subprocess_with_logging = subprocess_runner


def _sync_benchmark_dependencies():
    _benchmark.print = print
    _benchmark.query_sacct_for_jobs = globals()["query_sacct_for_jobs"]
    _benchmark.update_launch_metadata = globals()["update_launch_metadata"]


def _sync_run_log_dependencies():
    _sync_benchmark_dependencies()
    _run_logs.print = print
    _run_logs.section = section
    _run_logs.is_snakemake_locked = globals()["is_snakemake_locked"]
    _run_logs.generate_run_benchmark = globals()["generate_run_benchmark"]


def _sync_update_dependencies():
    _update_command.subprocess = subprocess
    _update_command.print = print
    version_func = globals()["get_installed_drakkar_version"]
    if _is_facade_function(version_func, "get_installed_drakkar_version"):
        version_func = _UPDATE_GET_INSTALLED_VERSION
    _update_command.get_installed_drakkar_version = version_func
    _update_command.display_update_success = globals()["display_update_success"]


class RichArgumentParser(_cli_help.RichArgumentParser):
    def print_help(self, file=None):
        _sync_help_dependencies()
        return super().print_help(file)

    def _print_message(self, message, file=None):
        _sync_help_dependencies()
        return super()._print_message(message, file=file)

    def error(self, message):
        _sync_help_dependencies()
        return super().error(message)


def _rich_available():
    _sync_help_dependencies()
    return _cli_help._rich_available()


def _rich_table(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_table(*args, **kwargs)


def _rich_action_table(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_action_table(*args, **kwargs)


def _rich_command_table(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_command_table(*args, **kwargs)


def _rich_subcommand_tables(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_subcommand_tables(*args, **kwargs)


def _rich_examples_panel(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_examples_panel(*args, **kwargs)


def _rich_action_tables(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_action_tables(*args, **kwargs)


def _rich_help_renderable(*args, **kwargs):
    _sync_help_dependencies()
    return _cli_help._rich_help_renderable(*args, **kwargs)


def normalize_annotation_type(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.normalize_annotation_type(*args, **kwargs)


def available_gtdb_versions(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.available_gtdb_versions(*args, **kwargs)


def validate_gtdb_version(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.validate_gtdb_version(*args, **kwargs)


def validate_database_version(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.validate_database_version(*args, **kwargs)


def validate_download_runtime(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.validate_download_runtime(*args, **kwargs)


def normalize_cataloging_binners(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.normalize_cataloging_binners(*args, **kwargs)


def positive_int(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.positive_int(*args, **kwargs)


def add_resource_multiplier_arguments(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.add_resource_multiplier_arguments(*args, **kwargs)


def add_benchmark_argument(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.add_benchmark_argument(*args, **kwargs)


def resource_config(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.resource_config(*args, **kwargs)


def default_resource_args(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.default_resource_args(*args, **kwargs)


def default_database_version(*args, **kwargs):
    _sync_validation_dependencies()
    return _VALIDATION_DEFAULT_DATABASE_VERSION(*args, **kwargs)


def validate_managed_database_version(*args, **kwargs):
    _sync_validation_dependencies()
    return _cli_validation.validate_managed_database_version(*args, **kwargs)


def replace_config_value(*args, **kwargs):
    _sync_config_dependencies()
    return _CONFIG_REPLACE_CONFIG_VALUE(*args, **kwargs)


def set_default_database_path(*args, **kwargs):
    _sync_config_dependencies()
    return _CONFIG_SET_DEFAULT_DATABASE_PATH(*args, **kwargs)


def resolve_editor_command(*args, **kwargs):
    return _CONFIG_RESOLVE_EDITOR_COMMAND(*args, **kwargs)


def view_config(*args, **kwargs):
    _sync_config_dependencies()
    return _CONFIG_VIEW_CONFIG(*args, **kwargs)


def edit_config(*args, **kwargs):
    _sync_config_dependencies()
    return _CONFIG_EDIT_CONFIG(*args, **kwargs)


def overwrite_output_directory(*args, **kwargs):
    return _OUTPUT_OVERWRITE_OUTPUT_DIRECTORY(*args, **kwargs)


def prompt_overwrite_locked_directory(*args, **kwargs):
    _sync_output_path_dependencies()
    return _OUTPUT_PROMPT_OVERWRITE_LOCKED_DIRECTORY(*args, **kwargs)


def prepare_output_directory(*args, **kwargs):
    _sync_output_path_dependencies()
    return _output_paths.prepare_output_directory(*args, **kwargs)


def validate_path(*args, **kwargs):
    _sync_output_path_dependencies()
    return _output_paths.validate_path(*args, **kwargs)


def validate_launch_metadata_directory(*args, **kwargs):
    _sync_output_path_dependencies()
    _run_metadata.validate_launch_metadata_directory = _output_paths.validate_launch_metadata_directory
    return _output_paths.validate_launch_metadata_directory(*args, **kwargs)


def write_launch_metadata(*args, **kwargs):
    _sync_output_path_dependencies()
    _run_metadata.validate_launch_metadata_directory = globals()["validate_launch_metadata_directory"]
    _run_metadata.print = print
    return _run_metadata.write_launch_metadata(*args, **kwargs)


def generate_run_benchmark(*args, **kwargs):
    _sync_benchmark_dependencies()
    return _benchmark.generate_run_benchmark(*args, **kwargs)


def run_logging(*args, **kwargs):
    _sync_run_log_dependencies()
    return _run_logs.run_logging(*args, **kwargs)


def run_unlock(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_unlock(*args, **kwargs)


def run_snakemake_environments(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_environments(*args, **kwargs)


def run_snakemake_preprocessing(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_preprocessing(*args, **kwargs)


def run_snakemake_cataloging(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_cataloging(*args, **kwargs)


def run_snakemake_cataloging2(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_cataloging2(*args, **kwargs)


def run_snakemake_profiling(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_profiling(*args, **kwargs)


def run_snakemake_dereplicating(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_dereplicating(*args, **kwargs)


def run_snakemake_annotating(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_annotating(*args, **kwargs)


def run_snakemake_inspecting(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_inspecting(*args, **kwargs)


def run_snakemake_expressing(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_expressing(*args, **kwargs)


def run_snakemake_database(*args, **kwargs):
    _sync_workflow_dependencies()
    return _workflow_launcher.run_snakemake_database(*args, **kwargs)


def get_installed_drakkar_version(*args, **kwargs):
    return _UPDATE_GET_INSTALLED_VERSION(*args, **kwargs)


def run_update(*args, **kwargs):
    _sync_update_dependencies()
    return _update_command.run_update(*args, **kwargs)


def _sync_main_dependencies():
    _cli_parser.RichArgumentParser = RichArgumentParser
    _cli_parser._set_help_metadata = _set_help_metadata
    _cli_main.build_parser = _cli_parser.build_parser
    _cli_main.config_vars = config_vars
    _cli_main.RichArgumentParser = RichArgumentParser
    _cli_main.prepare_output_directory = prepare_output_directory
    _cli_main.write_launch_metadata = write_launch_metadata
    _cli_main.run_logging = run_logging
    _cli_main.run_update = run_update
    _cli_main.run_sftp_transfer = run_sftp_transfer
    _cli_main.run_unlock = run_unlock
    _cli_main.run_snakemake_environments = run_snakemake_environments
    _cli_main.run_snakemake_preprocessing = run_snakemake_preprocessing
    _cli_main.run_snakemake_cataloging = run_snakemake_cataloging
    _cli_main.run_snakemake_profiling = run_snakemake_profiling
    _cli_main.run_snakemake_dereplicating = run_snakemake_dereplicating
    _cli_main.run_snakemake_annotating = run_snakemake_annotating
    _cli_main.run_snakemake_inspecting = run_snakemake_inspecting
    _cli_main.run_snakemake_expressing = run_snakemake_expressing
    _cli_main.run_snakemake_database = run_snakemake_database
    _cli_main.validate_path = validate_path
    _cli_main.validate_download_runtime = validate_download_runtime
    _cli_main.validate_gtdb_version = validate_gtdb_version
    _cli_main.validate_managed_database_version = validate_managed_database_version
    _cli_main.normalize_cataloging_binners = normalize_cataloging_binners
    _cli_main.normalize_annotation_type = normalize_annotation_type
    _cli_main.default_resource_args = default_resource_args
    _cli_main.resource_config = resource_config
    _cli_main.view_config = view_config
    _cli_main.edit_config = edit_config
    _cli_main.set_default_database_path = set_default_database_path


def main():
    _sync_main_dependencies()
    return _cli_main.main()


if __name__ == "__main__":
    main()
