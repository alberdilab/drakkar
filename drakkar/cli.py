import argparse
import os
import sys
import subprocess
import yaml
import re
import json
import shlex
import shutil
import pandas as pd
try:
    from importlib.metadata import PackageNotFoundError, version as get_distribution_version
except ImportError:  # pragma: no cover - Python < 3.8 fallback
    try:
        from importlib_metadata import PackageNotFoundError, version as get_distribution_version
    except ImportError:  # pragma: no cover - fallback if backport is absent
        PackageNotFoundError = Exception
        get_distribution_version = None
from pathlib import Path
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import PurePosixPath
try:
    from rich import box as rich_box
    from rich.console import Group as RichGroup
    from rich.panel import Panel as RichPanel
    from rich.table import Table as RichTable
    from rich.text import Text as RichText
except ImportError:  # pragma: no cover - only used if runtime deps are broken.
    rich_box = None
    RichGroup = None
    RichPanel = None
    RichTable = None
    RichText = None
from drakkar import __version__
from drakkar.database_registry import (
    MANAGED_DATABASES,
    database_release_dir,
    normalize_managed_database_name,
)
from drakkar.utils import *
from drakkar.output import get_console, print, prompt, section

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
# Define text colors
###

HEADER1 = "\033[1;95m"
ERROR = "\033[1;31m"
INFO = "\033[1;34m"
RESET = "\033[0m"

WORKFLOW_RUN_COMMANDS = {
    "complete",
    "preprocessing",
    "cataloging",
    "profiling",
    "dereplicating",
    "annotating",
    "inspecting",
    "expressing",
    "database",
    "environments",
}

READ_ONLY_COMMANDS = {"config", "logging"}

###
# Define helper functions
###

def normalize_annotation_type(annotation_type):
    functional_components = {
        "kegg", "cazy", "pfam", "virulence", "amr", "signalp",
        "dbcan", "antismash", "defense", "mobile"
    }
    gene_components = {"kegg", "cazy", "pfam", "virulence", "amr", "signalp"}
    aliases = {"vfdb": "virulence", "genomad": "mobile"}
    allowed = {
        "taxonomy", "function", "genes", "network",
        *functional_components
    }
    option_order = [
        "taxonomy", "function", "genes", "network",
        "kegg", "cazy", "pfam", "virulence", "amr", "signalp",
        "dbcan", "antismash", "defense", "mobile"
    ]
    items = [aliases.get(item.strip().lower(), item.strip().lower()) for item in annotation_type.split(",") if item.strip()]
    invalid = [item for item in items if item not in allowed]
    if not items or invalid:
        print(f"{ERROR}ERROR:{RESET} --annotation-type must be a comma-separated list including taxonomy, function, genes, kegg, cazy, pfam, virulence, amr, signalp, dbcan, antismash, defense, mobile, and/or network.")
        return None

    expanded = set(items)
    if "function" in expanded:
        expanded.update(functional_components)
    if "genes" in expanded:
        expanded.update(gene_components)

    normalized = [opt for opt in option_order if opt in expanded]
    return ",".join(normalized)


def available_gtdb_versions(config=None):
    source = config if config is not None else config_vars
    source = source or {}
    versions = []
    for key, value in source.items():
        match = re.fullmatch(r"GTDB_DB_(\d+)", str(key))
        if match and value:
            versions.append(match.group(1))
    return sorted(set(versions), key=lambda version: int(version), reverse=True)


def validate_gtdb_version(version, config=None):
    if not version:
        return None
    version = str(version).strip()
    if not re.fullmatch(r"\d+", version):
        print(f"{ERROR}ERROR:{RESET} --gtdb-version must be a GTDB release number such as 232.")
        return None

    versions = available_gtdb_versions(config)
    if version not in versions:
        supported = ", ".join(versions) if versions else "none configured"
        print(f"{ERROR}ERROR:{RESET} --gtdb-version {version} is not configured. Available versions: {supported}.")
        return None
    return version


def validate_database_version(version):
    version = (version or "").strip()
    if not version or version in {".", ".."} or "/" in version or "\\" in version:
        print(f"{ERROR}ERROR:{RESET} --version must be a single folder name, not a path.")
        return None
    return version


def validate_download_runtime(value):
    try:
        runtime = int(value)
    except (TypeError, ValueError):
        print(f"{ERROR}ERROR:{RESET} --download-runtime must be a positive integer number of minutes.")
        return None
    if runtime <= 0:
        print(f"{ERROR}ERROR:{RESET} --download-runtime must be a positive integer number of minutes.")
        return None
    return runtime


def positive_int(value):
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def add_resource_multiplier_arguments(parser):
    parser.add_argument(
        "--memory-multiplier",
        type=positive_int,
        default=1,
        help=(
            "Multiply Snakemake memory requests before applying SNAKEMAKE_MAX_GB. "
            "Default: 1"
        ),
    )
    parser.add_argument(
        "--time-multiplier",
        type=positive_int,
        default=1,
        help=(
            "Multiply Snakemake runtime requests before applying SNAKEMAKE_MAX_TIME. "
            "Default: 1"
        ),
    )


def resource_config(memory_multiplier=1, time_multiplier=1):
    return f"memory_multiplier={memory_multiplier} time_multiplier={time_multiplier} "


def default_resource_args(memory_multiplier=1, time_multiplier=1):
    max_mem_mb = int(config_vars.get("SNAKEMAKE_MAX_GB", 1024)) * 1024
    max_runtime = int(config_vars.get("SNAKEMAKE_MAX_TIME", 14 * 24 * 60))
    default_mem_mb = min(max_mem_mb, 8 * 1024 * memory_multiplier)
    default_runtime = min(max_runtime, 10 * time_multiplier)
    return f"--default-resources mem_mb={default_mem_mb} runtime={default_runtime} "


def default_database_version(database_name):
    if database_name == "vfdb":
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")
    return None


def validate_managed_database_version(database_name, version):
    if not version:
        default_version = default_database_version(database_name)
        if default_version:
            return default_version
        print(f"{ERROR}ERROR:{RESET} --version is required for {database_name}.")
        return None
    version = validate_database_version(version)
    if not version:
        return None
    if database_name == "kegg":
        try:
            parsed = datetime.strptime(version, "%Y-%m-%d")
        except ValueError:
            print(f"{ERROR}ERROR:{RESET} KEGG --version must be an archive date in YYYY-MM-DD format, e.g. 2026-02-01")
            return None
        return parsed.strftime("%Y-%m-%d")
    return version


def replace_config_value(config_key, new_value):
    pattern = re.compile(rf"^({re.escape(config_key)}:\s*)(\".*?\"|'.*?'|[^\n#]+)(\s*(#.*)?)?$", re.MULTILINE)
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    replacement = rf'\1"{new_value}"\3'
    updated_text, count = pattern.subn(replacement, config_text, count=1)
    if count != 1:
        raise ValueError(f"Could not update {config_key} in {CONFIG_PATH}")
    CONFIG_PATH.write_text(updated_text, encoding="utf-8")


def set_default_database_path(database_name, directory, version):
    definition = MANAGED_DATABASES[database_name]
    default_path = str(database_release_dir(database_name, directory, version))
    replace_config_value(definition["config_key"], default_path)
    return default_path


def resolve_editor_command():
    for env_var in ("VISUAL", "EDITOR"):
        value = os.environ.get(env_var)
        if value:
            return shlex.split(value)
    for candidate in ("nano", "vim", "vi"):
        resolved = shutil.which(candidate)
        if resolved:
            return [resolved]
    return None


def view_config():
    if not CONFIG_PATH.exists():
        print(f"{ERROR}ERROR:{RESET} config.yaml not found: {CONFIG_PATH}")
        return 1
    print(CONFIG_PATH.resolve())
    print("")
    text = CONFIG_PATH.read_text(encoding="utf-8")
    sys.stdout.write(text)
    if text and not text.endswith("\n"):
        sys.stdout.write("\n")
    return 0


def edit_config():
    if not CONFIG_PATH.exists():
        print(f"{ERROR}ERROR:{RESET} config.yaml not found: {CONFIG_PATH}")
        return 1
    editor_cmd = resolve_editor_command()
    if not editor_cmd:
        print(f"{ERROR}ERROR:{RESET} No terminal editor found. Set $VISUAL or $EDITOR.")
        return 1
    try:
        subprocess.run([*editor_cmd, str(CONFIG_PATH)], check=True)
    except FileNotFoundError:
        print(f"{ERROR}ERROR:{RESET} Editor not found: {' '.join(editor_cmd)}")
        return 1
    except subprocess.CalledProcessError as exc:
        print(f"{ERROR}ERROR:{RESET} Editor exited with code {exc.returncode}")
        return exc.returncode or 1
    return 0


def overwrite_output_directory(output_dir):
    output_path = Path(output_dir)
    if output_path.exists():
        shutil.rmtree(output_path)


def prompt_overwrite_locked_directory(output_dir):
    print(f"{INFO}INFO:{RESET} The output directory is locked: {output_dir}")
    print("This usually means a previous Snakemake run ended unexpectedly.")
    print(f"Use 'drakkar logging -o {output_dir}' to inspect the latest Snakemake log before unlocking or overwriting.")
    print("Overwriting will delete the entire directory and rerun from scratch.")
    while True:
        response = prompt("Delete the locked directory and overwrite it? [y/N]: ").strip().lower()
        if response in {"y", "yes"}:
            return True
        if response in {"", "n", "no"}:
            return False
        print("Please answer y or n.")


def _rich_available():
    return all((rich_box, RichGroup, RichPanel, RichTable, RichText))


def _clean_usage(parser):
    try:
        formatter = parser._get_formatter()
        formatter._width = 68
        formatter.add_usage(
            parser.usage,
            parser._actions,
            parser._mutually_exclusive_groups,
        )
        usage = formatter.format_help().strip()
    except Exception:
        usage = parser.format_usage().strip()
    if usage.startswith("usage: "):
        usage = usage[len("usage: "):]
    return usage


def _action_value(action):
    if getattr(action, "nargs", None) == 0:
        return ""
    if action.metavar is not None:
        if isinstance(action.metavar, tuple):
            value = " ".join(str(item) for item in action.metavar)
        else:
            value = str(action.metavar)
    elif action.choices is not None:
        value = "{" + ",".join(str(choice) for choice in action.choices) + "}"
    else:
        value = str(action.dest).upper()

    if action.nargs == "+":
        return f"{value} ..."
    if action.nargs == "*":
        return f"[{value} ...]"
    if action.nargs == "?":
        return f"[{value}]"
    return value


def _expanded_action_help(parser, action):
    if action.help in (None, argparse.SUPPRESS):
        return ""
    try:
        help_text = parser._get_formatter()._expand_help(action)
    except Exception:
        help_text = str(action.help)
    if getattr(action, "required", False):
        help_text = f"{help_text} Required." if help_text else "Required."
    return help_text


def _display_group_title(title):
    if title in {"options", "optional arguments"}:
        return "Options"
    if title == "positional arguments":
        return "Arguments"
    return str(title).replace("_", " ").title()


def _set_help_metadata(parser, *, category=None, examples=None, sections=None, command_groups=None):
    parser._help_category = category
    parser._help_examples = list(examples or [])
    parser._help_sections = list(sections or [])
    parser._help_command_groups = list(command_groups or [])
    return parser


def _rich_table(title):
    table = RichTable(
        title=RichText(title, style="drakkar.heading"),
        box=rich_box.SIMPLE_HEAVY,
        border_style="drakkar.rule",
        header_style="drakkar.heading",
        expand=True,
        padding=(0, 1),
    )
    table.add_column("Option", style="bold #d6a642", no_wrap=True)
    table.add_column("Value", style="#5f9ea0", no_wrap=True)
    table.add_column("Description", style="drakkar.help")
    return table


def _find_action_by_key(parser, key):
    for action in parser._actions:
        if action.dest == key:
            return action
        if key in getattr(action, "option_strings", ()):
            return action
    return None


def _rich_action_table(parser, action_group):
    actions = [
        action
        for action in action_group._group_actions
        if action.help is not argparse.SUPPRESS
        and not isinstance(action, argparse._SubParsersAction)
    ]
    if not actions:
        return None

    table = _rich_table(_display_group_title(action_group.title))
    for action in actions:
        option = ", ".join(action.option_strings) if action.option_strings else action.dest
        table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
    return table


def _subparser_action(parser):
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None


def _rich_command_table(title, rows):
    if not rows:
        return None
    table = RichTable(
        title=RichText(title, style="drakkar.heading"),
        box=rich_box.SIMPLE_HEAVY,
        border_style="drakkar.rule",
        header_style="drakkar.heading",
        expand=True,
        padding=(0, 1),
    )
    table.add_column("Command", style="bold #d6a642", no_wrap=True)
    table.add_column("Purpose", style="drakkar.help")
    for command, purpose in rows:
        table.add_row(command, purpose or "")
    return table


def _rich_subcommand_tables(parser):
    action = _subparser_action(parser)
    if action is None:
        return []

    rows_by_name = {}
    for subaction in action._get_subactions():
        command_name = str(subaction.metavar)
        rows_by_name[command_name] = (command_name, subaction.help or "")

    tables = []
    grouped_names = set()
    for title, names in getattr(parser, "_help_command_groups", []):
        rows = []
        for name in names:
            matched_name = None
            if name in rows_by_name:
                matched_name = name
            else:
                for row_name in rows_by_name:
                    if row_name.startswith(f"{name} "):
                        matched_name = row_name
                        break
            if matched_name is not None:
                rows.append(rows_by_name[matched_name])
                grouped_names.add(matched_name)
        table = _rich_command_table(title, rows)
        if table is not None:
            tables.append(table)

    remaining_rows = [rows_by_name[name] for name in rows_by_name if name not in grouped_names]
    remaining_table = _rich_command_table("Commands", remaining_rows)
    if remaining_table is not None:
        tables.append(remaining_table)
    return tables


def _rich_examples_panel(parser):
    examples = getattr(parser, "_help_examples", None) or []
    if not examples:
        return None

    body = RichText()
    for index, example in enumerate(examples):
        if index:
            body.append("\n")
        body.append("$ ", style="drakkar.heading")
        body.append(example, style="drakkar.help")
    return RichPanel(
        body,
        title=RichText(" Examples ", style="drakkar.heading"),
        border_style="drakkar.rule",
        padding=(1, 2),
    )


def _rich_action_tables(parser):
    custom_sections = getattr(parser, "_help_sections", None) or []
    if not custom_sections:
        return [
            table
            for action_group in parser._action_groups
            for table in [_rich_action_table(parser, action_group)]
            if table is not None
        ]

    tables = []
    assigned = set()
    for title, keys in custom_sections:
        actions = []
        for key in keys:
            action = _find_action_by_key(parser, key)
            if (
                action is None
                or action in assigned
                or action.help is argparse.SUPPRESS
                or isinstance(action, argparse._SubParsersAction)
            ):
                continue
            actions.append(action)
            assigned.add(action)
        if not actions:
            continue
        table = _rich_table(title)
        for action in actions:
            option = ", ".join(action.option_strings) if action.option_strings else action.dest
            table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
        tables.append(table)

    remaining = [
        action
        for action in parser._actions
        if action not in assigned
        and action.help is not argparse.SUPPRESS
        and not isinstance(action, argparse._SubParsersAction)
    ]
    if remaining:
        table = _rich_table("General")
        for action in remaining:
            option = ", ".join(action.option_strings) if action.option_strings else action.dest
            table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
        tables.append(table)
    return tables


def _rich_help_renderable(parser):
    if not _rich_available():
        return None

    header = RichText()
    if parser.description:
        header.append(parser.description, style="drakkar.text")
        header.append("\n\n")
    category = getattr(parser, "_help_category", None)
    if category:
        header.append("Category\n", style="drakkar.heading")
        header.append(str(category), style="drakkar.help")
        header.append("\n\n")
    header.append("Usage\n", style="drakkar.heading")
    header.append(_clean_usage(parser), style="drakkar.help")

    renderables = [
        RichPanel(
            header,
            title=RichText(f" {parser.prog} ", style="drakkar.heading"),
            border_style="drakkar.rule",
            padding=(1, 2),
        )
    ]

    renderables.extend(_rich_subcommand_tables(parser))

    examples_panel = _rich_examples_panel(parser)
    if examples_panel is not None:
        renderables.append(examples_panel)

    renderables.extend(_rich_action_tables(parser))

    return RichGroup(*renderables)


class RichArgumentParser(argparse.ArgumentParser):
    """Render argparse output through the shared Rich console."""

    def print_help(self, file=None):
        renderable = _rich_help_renderable(self)
        target = get_console(file)
        if target is not None and renderable is not None:
            target.print(renderable)
            return
        super().print_help(file)

    def _print_message(self, message, file=None):
        if message:
            print(message, end="", file=file, style="drakkar.help")

    def error(self, message):
        self.print_usage(sys.stderr)
        print(f"ERROR: {message}", file=sys.stderr)
        self.exit(2)


def prepare_output_directory(output_dir, overwrite=False):
    output_path = Path(output_dir)
    if not is_snakemake_locked(str(output_path)):
        return True

    if overwrite:
        print(f"{INFO}INFO:{RESET} Removing locked output directory and overwriting: {output_path}")
        overwrite_output_directory(output_path)
        return True

    if sys.stdin.isatty():
        if prompt_overwrite_locked_directory(output_path):
            overwrite_output_directory(output_path)
            return True
    else:
        print(f"{ERROR}ERROR:{RESET} Output directory is locked and no interactive prompt is available.")

    print(f"{ERROR}ERROR:{RESET} Output directory remains locked: {output_path}")
    print(f"{INFO}Use 'drakkar logging -o {output_path}' to inspect the latest Snakemake log, --overwrite to delete it automatically, or 'drakkar unlock -o {output_path}' if you only want to clear the Snakemake lock.{RESET}")
    return False

def validate_path(path_value, label, expect_dir=False, allow_url=False):
    if not path_value:
        return True
    if allow_url and not expect_dir and is_url(path_value):
        return True
    if expect_dir:
        exists = os.path.isdir(path_value)
    else:
        exists = os.path.exists(path_value)
    if not exists:
        kind = "directory" if expect_dir else "path"
        print(f"{ERROR}ERROR:{RESET} {label} {kind} not found: {path_value}")
        return False
    return True

def validate_launch_metadata_directory(output_dir):
    output_path = Path(output_dir)
    existing_path = output_path
    while not existing_path.exists() and existing_path.parent != existing_path:
        existing_path = existing_path.parent

    if not existing_path.exists():
        print(f"{ERROR}ERROR:{RESET} Cannot find a parent directory for output path: {output_path}")
        return False
    if not existing_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path parent is not a directory: {existing_path}")
        return False

    access_target = output_path if output_path.exists() else existing_path
    if output_path.exists() and not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path is not a directory: {output_path}")
        return False
    if not os.access(access_target, os.R_OK | os.W_OK | os.X_OK):
        print(f"{ERROR}ERROR:{RESET} Cannot write Drakkar run metadata in: {output_path}")
        print("The current/output directory must be readable, writable, and searchable.")
        print("Run drakkar from a writable directory or pass -o/--output to a writable output directory.")
        return False

    return True

def get_modules_to_run(command):
    if command == "complete":
        return ["preprocessing", "cataloging", "profiling", "annotating"]
    if command:
        return [command]
    return []

def build_snakemake_log_path(output_dir, run_id):
    return Path(output_dir) / "log" / f"drakkar_{run_id}.snakemake.log"


def write_launch_metadata(args, output_dir, env_path=None):
    output_path = Path(output_dir)
    if not validate_launch_metadata_directory(output_path):
        return None
    try:
        output_path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"{ERROR}ERROR:{RESET} Cannot create output directory for Drakkar run metadata: {output_path}")
        print(f"{exc.__class__.__name__}: {exc}")
        return None
    timestamp = datetime.now(timezone.utc)
    run_id = timestamp.strftime("%Y%m%d-%H%M%S")
    snakemake_log_path = None
    if args.command in WORKFLOW_RUN_COMMANDS:
        snakemake_log_path = build_snakemake_log_path(output_path, run_id)
        try:
            snakemake_log_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"{ERROR}ERROR:{RESET} Cannot create log directory for Drakkar run metadata: {snakemake_log_path.parent}")
            print(f"{exc.__class__.__name__}: {exc}")
            return None
    metadata = {
        "run_id": run_id,
        "timestamp": timestamp.isoformat(),
        "started_at": timestamp.isoformat(),
        "command": args.command,
        "modules": get_modules_to_run(args.command),
        "working_directory": str(Path.cwd()),
        "output_directory": str(output_path.resolve()),
        "arguments": vars(args),
        "argv": sys.argv,
        "status": "prepared",
    }
    if env_path is not None:
        metadata["env_path"] = env_path
    if snakemake_log_path is not None:
        metadata["snakemake_log"] = str(snakemake_log_path.resolve())
    metadata_path = output_path / f"drakkar_{run_id}.yaml"
    try:
        with open(metadata_path, "w") as f:
            yaml.safe_dump(metadata, f, sort_keys=False)
    except OSError as exc:
        print(f"{ERROR}ERROR:{RESET} Cannot write Drakkar run metadata: {metadata_path}")
        print("Run drakkar from a writable directory or pass -o/--output to a writable output directory.")
        print(f"{exc.__class__.__name__}: {exc}")
        return None
    return {
        "run_id": run_id,
        "metadata_path": metadata_path,
        "snakemake_log_path": snakemake_log_path,
    }


def update_launch_metadata(metadata_path, **updates):
    if not metadata_path:
        return None
    metadata_path = Path(metadata_path)
    if not metadata_path.exists():
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as handle:
            metadata = yaml.safe_load(handle) or {}
    except OSError:
        return None
    metadata.update(updates)
    try:
        with open(metadata_path, "w", encoding="utf-8") as handle:
            yaml.safe_dump(metadata, handle, sort_keys=False)
    except OSError:
        return None
    return metadata


def finalize_launch_metadata(run_info, status, exit_code=None, current_workflow=None):
    if not run_info:
        return None
    payload = {
        "status": status,
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    if exit_code is not None:
        payload["exit_code"] = exit_code
    if current_workflow is not None:
        payload["current_workflow"] = current_workflow
    return update_launch_metadata(run_info["metadata_path"], **payload)


def run_subprocess_with_logging(command, run_info=None, workflow_name=None):
    metadata_path = run_info["metadata_path"] if run_info else None
    log_path = Path(run_info["snakemake_log_path"]) if run_info and run_info.get("snakemake_log_path") else None
    if metadata_path:
        update_launch_metadata(
            metadata_path,
            status="running",
            current_workflow=workflow_name,
        )

    log_handle = None
    try:
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_handle = open(log_path, "a", encoding="utf-8")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        if process.stdout is not None:
            try:
                for line in process.stdout:
                    print(line, end="")
                    if log_handle is not None:
                        log_handle.write(line)
            finally:
                process.stdout.close()
        return_code = process.wait()
    except Exception:
        if log_handle is not None:
            log_handle.flush()
            log_handle.close()
        finalize_launch_metadata(run_info, "failed", current_workflow=workflow_name)
        raise
    finally:
        if log_handle is not None and not log_handle.closed:
            log_handle.flush()
            log_handle.close()

    if return_code != 0:
        finalize_launch_metadata(run_info, "failed", return_code, current_workflow=workflow_name)
        raise subprocess.CalledProcessError(return_code, command)

    finalize_launch_metadata(run_info, "success", 0, current_workflow=workflow_name)
    return return_code


def workflow_run_sort_key(metadata_path):
    path = Path(metadata_path)
    match = re.search(r"drakkar_(\d{8}-\d{6})\.ya?ml$", path.name)
    if match:
        return match.group(1)
    return path.name


def discover_run_metadata(output_dir):
    output_path = Path(output_dir)
    runs = []
    for metadata_path in sorted(output_path.glob("drakkar_*.yaml"), key=workflow_run_sort_key, reverse=True):
        try:
            with open(metadata_path, "r", encoding="utf-8") as handle:
                metadata = yaml.safe_load(handle) or {}
        except OSError:
            continue
        command = metadata.get("command")
        if command not in WORKFLOW_RUN_COMMANDS:
            continue
        runs.append((metadata_path, metadata))
    return runs


def resolve_run_metadata(output_dir, run_id=None):
    runs = discover_run_metadata(output_dir)
    if not runs:
        return None, None
    if not run_id:
        return runs[0]

    selector = str(run_id).strip()
    normalized_selector = selector.removeprefix("drakkar_").removesuffix(".yaml")
    for metadata_path, metadata in runs:
        run_value = str(metadata.get("run_id", "")).strip()
        if selector == metadata_path.name or normalized_selector == run_value:
            return metadata_path, metadata
    return None, None


def discover_snakemake_fallback_logs(output_dir):
    log_dir = Path(output_dir) / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    return sorted((path for path in log_dir.glob("*") if path.is_file()), key=lambda path: path.stat().st_mtime, reverse=True)


def tail_file(path, line_count):
    lines = deque(maxlen=line_count)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lines.append(line.rstrip("\n"))
    return list(lines)


def extract_failure_excerpt(path, line_count=40):
    markers = (
        "RuleException",
        "MissingInputException",
        "WorkflowError",
        "CalledProcessError",
        "LockException",
        "Error in rule",
        "Traceback (most recent call last):",
    )
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle]
    last_index = None
    for index, line in enumerate(lines):
        if any(marker in line for marker in markers):
            last_index = index
    if last_index is None:
        return []
    start = max(0, last_index - 5)
    end = min(len(lines), last_index + line_count)
    return lines[start:end]


def classify_error_line(line):
    known_error_types = (
        "RuleException",
        "MissingInputException",
        "WorkflowError",
        "CalledProcessError",
        "LockException",
        "InputFunctionException",
        "ChildIOException",
    )
    if line.startswith("Error in rule "):
        return "RuleError"
    for error_type in known_error_types:
        if re.search(rf"\b{re.escape(error_type)}\b", line):
            return error_type
    match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*(?:Exception|Error))(?::|\b)", line)
    if match:
        return match.group(1)
    return None


def summarize_snakemake_log(path, metadata=None):
    summary = {
        "planned_jobs": None,
        "completed_steps": None,
        "total_steps": None,
        "percent_complete": None,
        "unique_rules": 0,
        "rule_executions": 0,
        "failed_rules": 0,
        "error_types": Counter(),
        "most_active_rules": [],
    }
    if path is None or not Path(path).exists():
        return summary

    rule_counter = Counter()
    finished_job_ids = set()
    started_job_ids = set()
    planned_jobs = None
    best_progress = None
    failed_rules = 0
    error_types = Counter()

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            match = re.match(r"^(?:local)?rule\s+(.+?):\s*$", stripped)
            if match:
                rule_counter[match.group(1)] += 1

            match = re.search(r"\bjobid:\s*(\d+)\b", stripped)
            if match:
                started_job_ids.add(match.group(1))

            match = re.search(r"Finished jobid:\s*(\d+)", stripped)
            if match:
                finished_job_ids.add(match.group(1))

            match = re.match(r"(\d+)\s+of\s+(\d+)\s+steps\s+\((\d+)%\)\s+done", stripped)
            if match:
                progress = (int(match.group(1)), int(match.group(2)), int(match.group(3)))
                if best_progress is None or progress[2] > best_progress[2] or (
                    progress[2] == best_progress[2] and progress[0] > best_progress[0]
                ):
                    best_progress = progress

            match = re.match(r"total\s+(\d+)\s*$", stripped)
            if match:
                count = int(match.group(1))
                planned_jobs = count if planned_jobs is None else max(planned_jobs, count)

            if stripped.startswith("Error in rule "):
                failed_rules += 1

            error_type = classify_error_line(stripped)
            if error_type:
                error_types[error_type] += 1

    completed_steps = best_progress[0] if best_progress else None
    total_steps = best_progress[1] if best_progress else None
    percent_complete = best_progress[2] if best_progress else None

    if completed_steps is None and finished_job_ids:
        completed_steps = len(finished_job_ids)
    if total_steps is None and planned_jobs is not None:
        total_steps = planned_jobs
    if percent_complete is None and completed_steps is not None and total_steps:
        percent_complete = int(round((completed_steps / total_steps) * 100))

    status = str((metadata or {}).get("status", "")).strip().lower()
    if status == "success":
        if total_steps is None:
            total_steps = planned_jobs or len(finished_job_ids) or len(started_job_ids) or None
        if total_steps is not None:
            completed_steps = total_steps
            percent_complete = 100

    summary.update(
        {
            "planned_jobs": planned_jobs,
            "completed_steps": completed_steps,
            "total_steps": total_steps,
            "percent_complete": percent_complete,
            "unique_rules": len(rule_counter),
            "rule_executions": sum(rule_counter.values()),
            "failed_rules": failed_rules,
            "error_types": error_types,
            "most_active_rules": rule_counter.most_common(5),
        }
    )
    return summary


def print_snakemake_summary(summary):
    section("EXECUTION SUMMARY")
    planned_jobs = summary.get("planned_jobs")
    completed_steps = summary.get("completed_steps")
    total_steps = summary.get("total_steps")
    percent_complete = summary.get("percent_complete")
    unique_rules = summary.get("unique_rules", 0)
    rule_executions = summary.get("rule_executions", 0)
    failed_rules = summary.get("failed_rules", 0)
    error_types = summary.get("error_types") or Counter()
    most_active_rules = summary.get("most_active_rules") or []

    print(f"Planned jobs: {planned_jobs if planned_jobs is not None else 'unknown'}")
    if percent_complete is not None and completed_steps is not None and total_steps is not None:
        print(f"Workflow progress: {percent_complete}% ({completed_steps}/{total_steps} steps)")
    elif completed_steps is not None:
        print(f"Completed steps observed: {completed_steps}")
    else:
        print("Workflow progress: unknown")
    print(f"Rules observed: {unique_rules} unique, {rule_executions} executions")
    print(f"Failed rules detected: {failed_rules}")
    if error_types:
        formatted_types = ", ".join(f"{name} ({count})" for name, count in error_types.most_common())
        print(f"Error types: {formatted_types}")
    else:
        print("Error types: none detected")
    if most_active_rules:
        formatted_rules = ", ".join(f"{name} ({count})" for name, count in most_active_rules)
        print(f"Most active rules: {formatted_rules}")


def run_logging(output_dir, run_id=None, tail=50, full=False, paths=False, list_runs=False, summary=False):
    output_path = Path(output_dir).resolve()
    if not output_path.exists():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output path is not a directory: {output_path}")
        return 1

    section("DRAKKAR LOGGING")
    print(f"Output directory: {output_path}")
    print(f"Locked: {'yes' if is_snakemake_locked(str(output_path)) else 'no'}")

    runs = discover_run_metadata(output_path)
    if list_runs:
        if not runs:
            print(f"{INFO}INFO:{RESET} No workflow run metadata found in {output_path}.")
        else:
            section("AVAILABLE RUNS")
            for metadata_path, metadata in runs:
                run_value = metadata.get("run_id", metadata_path.stem.removeprefix("drakkar_"))
                command = metadata.get("command", "unknown")
                status = metadata.get("status", "unknown")
                print(f"{run_value}: command={command}, status={status}")
        return 0

    metadata_path, metadata = resolve_run_metadata(output_path, run_id)
    if run_id and metadata is None:
        print(f"{ERROR}ERROR:{RESET} Run not found in {output_path}: {run_id}")
        return 1
    snakemake_log_path = None
    fallback_logs = discover_snakemake_fallback_logs(output_path)
    if metadata is not None:
        configured_log = metadata.get("snakemake_log")
        if configured_log:
            snakemake_log_path = Path(configured_log)
        elif metadata.get("run_id"):
            candidate = build_snakemake_log_path(output_path, metadata["run_id"])
            if candidate.exists():
                snakemake_log_path = candidate
    if snakemake_log_path is None and fallback_logs:
        snakemake_log_path = fallback_logs[0]

    if metadata is not None:
        section("RUN SUMMARY")
        print(f"Run ID: {metadata.get('run_id', metadata_path.stem.removeprefix('drakkar_'))}")
        print(f"Command: {metadata.get('command', 'unknown')}")
        modules = metadata.get("modules") or []
        print(f"Modules: {', '.join(modules) if modules else 'unknown'}")
        print(f"Status: {metadata.get('status', 'unknown')}")
        print(f"Started: {metadata.get('started_at', metadata.get('timestamp', 'unknown'))}")
        if metadata.get("finished_at"):
            print(f"Finished: {metadata['finished_at']}")
        if "exit_code" in metadata:
            print(f"Exit code: {metadata['exit_code']}")
        if metadata.get("current_workflow"):
            print(f"Current workflow: {metadata['current_workflow']}")
        print(f"Metadata file: {metadata_path}")
    else:
        print(f"{INFO}INFO:{RESET} No workflow metadata found. Falling back to Snakemake logs only.")

    log_summary = summarize_snakemake_log(snakemake_log_path, metadata=metadata)
    if snakemake_log_path is not None and snakemake_log_path.exists():
        print_snakemake_summary(log_summary)

    if paths:
        section("LOG PATHS")
        if snakemake_log_path is not None:
            print(f"Main Snakemake log: {snakemake_log_path}")
        else:
            print("Main Snakemake log: not found")
        for fallback_log in fallback_logs[:5]:
            print(f"Snakemake fallback log: {fallback_log}")
        extra_logs = sorted(
            path for path in (output_path / "log").rglob("*")
            if path.is_file() and path != snakemake_log_path
        ) if (output_path / "log").exists() else []
        for extra_log in extra_logs[:20]:
            print(f"Additional log: {extra_log}")

    if snakemake_log_path is None or not snakemake_log_path.exists():
        print(f"{INFO}INFO:{RESET} No Snakemake log file found in {output_path}.")
        return 0

    if summary and not full:
        return 0

    section("SNAKEMAKE LOG")
    print(f"Log file: {snakemake_log_path}")
    if full:
        with open(snakemake_log_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                print(line, end="")
        return 0

    excerpt = extract_failure_excerpt(snakemake_log_path)
    if excerpt:
        print("Most recent failure excerpt:")
        for line in excerpt:
            print(line)
        return 0

    print(f"Last {tail} lines:")
    for line in tail_file(snakemake_log_path, tail):
        print(line)
    return 0

def normalize_genome_name(name):
    if not name:
        return ""
    base = os.path.basename(str(name).strip())
    if base.endswith(".gz"):
        base = base[:-3]
    for ext in (".fa", ".fna", ".fasta"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base

def load_bins_map(output_dir):
    bins_path = Path(output_dir) / "data" / "bins_to_files.json"
    if not bins_path.exists():
        return {}
    with open(bins_path, "r") as f:
        return json.load(f)

def validate_and_write_quality_file(quality_path, output_dir):
    if not quality_path:
        return False
    if not os.path.isfile(quality_path):
        print(f"{ERROR}ERROR:{RESET} Quality file not found: {quality_path}")
        return False

    df = pd.read_csv(quality_path, sep=None, engine="python", encoding="utf-8-sig")
    col_map = {c: str(c).strip().lstrip("\ufeff").lower() for c in df.columns}
    df.rename(columns=col_map, inplace=True)
    required = {"genome", "completeness", "contamination"}
    if not required.issubset(set(df.columns)):
        print(f"{ERROR}ERROR:{RESET} Quality file must contain columns: genome, completeness, contamination")
        return False

    bins_map = load_bins_map(output_dir)
    if not bins_map:
        print(f"{ERROR}ERROR:{RESET} bins_to_files.json not found; cannot validate quality file.")
        return False

    # Build mapping from stem and basename to basename-with-extension
    stem_to_base = {}
    base_set = set()
    for path in bins_map.values():
        base = os.path.basename(path)
        if base.endswith(".gz"):
            base = base[:-3]
        base_set.add(base)
        stem_to_base[normalize_genome_name(base)] = base

    mapped = []
    missing = []
    for name in df["genome"].astype(str):
        key = name.strip()
        if key in base_set:
            mapped.append(key)
            continue
        stem = normalize_genome_name(key)
        mapped_base = stem_to_base.get(stem)
        if mapped_base:
            mapped.append(mapped_base)
        else:
            mapped.append(key)
            missing.append(key)

    if missing:
        print(f"{ERROR}ERROR:{RESET} Quality file missing bins: {missing}")
        return False

    out_dir = Path(output_dir) / "cataloging" / "final"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "all_bin_metadata.csv"
    out_df = df.copy()
    out_df["genome"] = mapped
    out_df = out_df[["genome", "completeness", "contamination"]]
    out_df.to_csv(out_path, index=False)
    return True
def list_files_recursive(base_dir, exclude_snakemake=False):
    for root, dirs, files in os.walk(base_dir):
        if exclude_snakemake:
            dirs[:] = [d for d in dirs if d != ".snakemake"]
        for filename in files:
            yield Path(root) / filename

def collect_transfer_files(base_dir, args):
    selected_files = set()
    missing_paths = []
    warnings = []

    def add_annotations():
        add_file(Path("annotating/cluster_annotations.tsv.xz"))
        add_file(Path("annotating/gene_annotations.tsv.xz"))
        add_file(Path("annotating/genome_taxonomy.tsv"))

    def add_mags():
        add_dir(Path("profiling_genomes/drep/dereplicated_genomes"))

    def add_profile():
        add_file(Path("profiling_genomes/final/bases.tsv"))
        add_file(Path("profiling_genomes/final/counts.tsv"))
        add_file(Path("profiling_genomes.tsv"))

    def add_expression():
        add_file(Path("expressing/gene_counts.tsv.xz"))

    def add_bins():
        add_dir(Path("cataloging/final"))

    def add_runner_logs():
        log_files = sorted(base_dir.glob("drakkar_*.yaml"))
        if log_files:
            selected_files.update(log_files)
        else:
            warnings.append("No drakkar_*.yaml run logs found in the local directory.")

    def add_file(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_file():
            selected_files.add(local_path)
        else:
            missing_paths.append(str(local_path))

    def add_dir(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_dir():
            for file_path in local_path.rglob("*"):
                if file_path.is_file():
                    selected_files.add(file_path)
        else:
            missing_paths.append(str(local_path))

    if args.all:
        for file_path in list_files_recursive(base_dir):
            if file_path.is_file():
                selected_files.add(file_path)
        add_runner_logs()
        return selected_files, missing_paths, warnings

    if args.data:
        for file_path in list_files_recursive(base_dir, exclude_snakemake=True):
            if file_path.is_file():
                selected_files.add(file_path)

    if args.results:
        add_annotations()
        add_mags()
        add_profile()
        add_expression()
        add_bins()

    if args.annotations:
        add_annotations()

    if args.mags:
        add_mags()

    if args.profile:
        add_profile()

    if args.expression:
        add_expression()

    if args.bins:
        add_bins()

    add_runner_logs()

    return selected_files, missing_paths, warnings

def build_sftp_batch_commands(files, base_dir, remote_dir):
    remote_root = PurePosixPath(remote_dir)
    mkdirs = set()
    put_cmds = []

    for file_path in sorted(files):
        rel_path = file_path.relative_to(base_dir)
        remote_path = remote_root / PurePosixPath(rel_path.as_posix())
        parent = remote_path.parent
        while parent != remote_root and parent not in mkdirs:
            mkdirs.add(parent)
            parent = parent.parent
        put_cmds.append((file_path, remote_path))

    commands = []
    for directory in sorted(mkdirs, key=lambda p: len(p.parts)):
        if directory == remote_root:
            continue
        commands.append(f'-mkdir "{directory.as_posix()}"')
    for local_path, remote_path in put_cmds:
        commands.append(f'put "{local_path}" "{remote_path.as_posix()}"')

    return "\n".join(commands) + "\n", put_cmds

def run_sftp_transfer(args):
    base_dir = Path(args.local_dir).resolve()
    if args.erda:
        args.host = "io.erda.dk"
        args.user = "antton.alberdi@snm.ku.dk"
    if not args.host or not args.user:
        print(f"{ERROR}ERROR:{RESET} --host and --user are required unless --erda is set.")
        return
    sftp_cmd = ["sftp"]
    if args.port:
        sftp_cmd += ["-P", str(args.port)]
    if args.identity:
        sftp_cmd += ["-i", args.identity]
    sftp_cmd += ["-b", "-", f"{args.user}@{args.host}"]
    check_cmds = f'ls "{args.remote_dir}"\n'
    check_result = subprocess.run(sftp_cmd, input=check_cmds, text=True, capture_output=True)
    if check_result.returncode != 0:
        stderr = check_result.stderr.strip()
        if stderr:
            print(stderr, file=sys.stderr)
        print(f"{ERROR}ERROR:{RESET} Remote directory not found: {args.remote_dir}")
        return
    selected_files, missing_paths, warnings = collect_transfer_files(base_dir, args)

    if warnings:
        for warning in warnings:
            print(f"{INFO}INFO:{RESET} {warning}")

    if missing_paths:
        for missing in missing_paths:
            print(f"{ERROR}MISSING:{RESET} {missing}")

    if not selected_files:
        print(f"{ERROR}ERROR:{RESET} No files selected for transfer.")
        return

    batch_commands, put_cmds = build_sftp_batch_commands(selected_files, base_dir, args.remote_dir)
    if args.verbose:
        for local_path, remote_path in put_cmds:
            print(f"{INFO}PUT:{RESET} {local_path} -> {remote_path.as_posix()}")

    result = subprocess.run(sftp_cmd, input=batch_commands, text=True, capture_output=True)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        mkdir_only_errors = False
        if stderr:
            error_lines = [line for line in stderr.splitlines() if line.strip()]
            mkdir_only_errors = all("Couldn't create directory" in line for line in error_lines)
        if mkdir_only_errors:
            print(f"{INFO}INFO:{RESET} sftp reported existing directories; continuing.")
        else:
            if stderr:
                print(stderr, file=sys.stderr)
            print(f"{ERROR}ERROR:{RESET} sftp transfer failed with exit code {result.returncode}")

###
# Define workflow launching functions
###

def run_unlock(workflow, output_dir, profile):

    unlock_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked")
    print(f"You can now rerun a new workflow using any drakkar command:")

def run_snakemake_environments(workflow, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    cmd = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {Path.cwd()} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} workflow={workflow} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(cmd, run_info=run_info, workflow_name=workflow)

def run_snakemake_preprocessing(
    workflow,
    project_name,
    output_dir,
    reference,
    env_path,
    profile,
    fraction=False,
    nonpareil=False,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):

    """ Run the preprocessing workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} reference={reference} fraction={fraction} nonpareil={nonpareil} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
        f"--slurm-delete-logfiles-older-than 0"
    ]

    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)


def run_snakemake_cataloging(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):

    """ Run the cataloging workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
    ]

    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

#Screen output control
def run_snakemake_cataloging2(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1):

    """ Run the cataloging workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--conda-frontend mamba "
        f"--use-conda "
    ]

    process = subprocess.Popen(snakemake_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    for line in process.stdout:
        print(line, end="")  # Print each line as it comes

    # Wait for the process to finish
    process.wait()

    if process.returncode != 0:
        # Read stderr and print it
        error_message = process.stderr.read()
        if "LockException" in error_message:
            display_unlock()
        else:
            print(f"\nERROR: Snakemake failed with exit code {process.returncode}!", file=sys.stderr)
            sys.exit(1)
    else:
        display_end()

def run_snakemake_profiling(
    workflow,
    project_name,
    profiling_type,
    output_dir,
    env_path,
    profile,
    fraction,
    ani,
    ignore_quality,
    quality_file,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the profiling workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} profiling_type={profiling_type} output_dir={output_dir} fraction={fraction} DREP_ANI={ani} "
        f"IGNORE_QUALITY={ignore_quality} QUALITY_FILE={bool(quality_file)} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_dereplicating(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    ani,
    ignore_quality,
    quality_file,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the dereplicating workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} DREP_ANI={ani} IGNORE_QUALITY={ignore_quality} QUALITY_FILE={bool(quality_file)} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_annotating(
    workflow,
    project_name,
    annotating_type,
    output_dir,
    env_path,
    profile,
    gtdb_version=None,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """ Run the profiling workflow """

    gtdb_config = f"gtdb_version={gtdb_version} " if gtdb_version else ""
    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} annotating_type={annotating_type} output_dir={output_dir} {gtdb_config}{resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_inspecting(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    """ Run the profiling workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c", 
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_expressing(workflow, project_name, output_dir, env_path, profile, memory_multiplier=1, time_multiplier=1, run_info=None):
    """ Run the expressing workflow """

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)

def run_snakemake_database(
    workflow,
    project_name,
    output_dir,
    env_path,
    profile,
    database_name,
    database_directory,
    database_version,
    download_runtime,
    memory_multiplier=1,
    time_multiplier=1,
    run_info=None,
):
    """Run a single database preparation workflow."""

    resource_overrides = resource_config(memory_multiplier, time_multiplier)
    default_resources = default_resource_args(memory_multiplier, time_multiplier)
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} project_name={project_name} workflow={workflow} output_dir={output_dir} "
        f"database_name={database_name} database_directory={database_directory} database_version={database_version} "
        f"database_download_runtime={download_runtime} {resource_overrides}"
        f"{default_resources}"
        f"--conda-prefix {env_path} "
        f"--use-conda "
    ]
    run_subprocess_with_logging(snakemake_command, run_info=run_info, workflow_name=workflow)


def get_installed_drakkar_version():
    if get_distribution_version is None:
        return __version__
    try:
        return get_distribution_version("drakkar")
    except PackageNotFoundError:
        return __version__
    except Exception:
        return __version__


def run_update():
    pip_cmd = [
        sys.executable, "-m", "pip", "install",
        "--upgrade", "--force-reinstall", "--no-deps",
        "git+https://github.com/alberdilab/drakkar.git",
    ]
    try:
        update_result = subprocess.run(pip_cmd)
    except Exception as exc:
        print(f"Update failed: {exc}", file=sys.stderr, flush=True)
        return 1
    if update_result.returncode != 0:
        return update_result.returncode
    display_update_success(get_installed_drakkar_version())
    return 0

###
# Main function to launch workflows
###

def main():
    parser = RichArgumentParser(
        prog="drakkar",
        description="Drakkar: A Snakemake-based workflow for sequencing analysis",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--version", action="version", version=f"drakkar {__version__}")
    subparsers = parser.add_subparsers(dest="command", help="Available workflows")

    # Define subcommands for each workflow
    subparser_complete = subparsers.add_parser("complete", help="Run the complete workflow")
    subparser_complete.add_argument("-i", "--input", required=False, help="Input directory")
    subparser_complete.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_complete.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    complete_reference_group = subparser_complete.add_mutually_exclusive_group()
    complete_reference_group.add_argument("-r", "--reference", required=False, help="Reference host genome FASTA")
    complete_reference_group.add_argument("-x", "--reference-index", required=False, help="Tarball containing a reference FASTA and Bowtie2 index files")
    subparser_complete.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
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
    add_resource_multiplier_arguments(subparser_complete)

    subparser_preprocessing = subparsers.add_parser("preprocessing", help="Run the preprocessing workflow (quality-filtering and host removal)")
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
    add_resource_multiplier_arguments(subparser_preprocessing)

    subparser_cataloging = subparsers.add_parser("cataloging", help="Run the cataloging (assembly and binning) workflow")
    subparser_cataloging.add_argument("-i", "--input", required=False, help="Input directory (required if no sample detail file is provided)")
    subparser_cataloging.add_argument("-f", "--file", required=False, help="Sample detail file (required if no input directory is provided)")
    subparser_cataloging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_cataloging.add_argument("-m", "--mode", required=False, help="Comma-separated list of cataloging modes (e.g. individual,all)")
    subparser_cataloging.add_argument("-c", "--multicoverage", action="store_true", help="Map samples sharing the same coverage group to each other's individual assemblies")
    subparser_cataloging.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_cataloging.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_cataloging.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_resource_multiplier_arguments(subparser_cataloging)

    subparser_profiling = subparsers.add_parser("profiling", help="Run the profiling workflow")
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
    add_resource_multiplier_arguments(subparser_profiling)

    subparser_dereplicating = subparsers.add_parser("dereplicating", help="Run dereplication only (no mapping)")
    subparser_dereplicating.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_dereplicating.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_dereplicating.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_dereplicating.add_argument("-a", "--ani", required=False, type=float, default=0.98, help="ANI threshold for dRep dereplication (-sa). Default: 0.98")
    subparser_dereplicating.add_argument("-n", "--ignore_quality", action="store_true", help="Pass --ignoreGenomeQuality to dRep during dereplication")
    subparser_dereplicating.add_argument("-q", "--quality", type=str, help="CSV/TSV with genome, completeness, contamination columns")
    subparser_dereplicating.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_dereplicating.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_dereplicating.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_resource_multiplier_arguments(subparser_dereplicating)

    subparser_annotating = subparsers.add_parser("annotating", help="Run the annotating workflow")
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
    add_resource_multiplier_arguments(subparser_annotating)

    subparser_inspecting = subparsers.add_parser("inspecting", help="Run the inspecting workflow")
    subparser_inspecting.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_inspecting.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_inspecting.add_argument("-m", "--mapping_dir", required=False, help="Directory containing the mapping (.bam) files")
    subparser_inspecting.add_argument("-c", "--cov_file", required=False, help="Tab-separated file containing mapping coverage per genome per smaple")
    subparser_inspecting.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_inspecting.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_inspecting.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_inspecting.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_resource_multiplier_arguments(subparser_inspecting)

    subparser_expressing = subparsers.add_parser("expressing", help="Run the microbial gene expression workflow")
    subparser_expressing.add_argument("-b", "--bins_dir", required=False, help="Directory in which bins (.fa or .fna) are stored")
    subparser_expressing.add_argument("-B", "--bins_file", required=False, help="Text file containing paths to the bins (.fa or .fna)")
    subparser_expressing.add_argument("-r", "--reads_dir", required=False, help="Directory in which metagenomic reads are stored")
    subparser_expressing.add_argument("-R", "--reads_file", required=False, help="Sample detail file")
    subparser_expressing.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_expressing.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_expressing.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    subparser_expressing.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    add_resource_multiplier_arguments(subparser_expressing)

    database_parent = RichArgumentParser(add_help=False)
    database_parent.add_argument("--directory", required=True, help="Base directory where the database release directory will be created")
    database_parent.add_argument("--version", required=False, help="Release folder name to create inside --directory (optional for vfdb; defaults to the UTC download date)")
    database_parent.add_argument("--download-runtime", type=int, default=120, help="Runtime in minutes for the database download/preparation rule (default: 120)")
    database_parent.add_argument("--overwrite", action="store_true", help="Delete a locked output directory and rerun from scratch")
    database_parent.add_argument("--set-default", action="store_true", help="Update config.yaml to use this installed database release by default")
    database_parent.add_argument("-e", "--env_path", type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    database_parent.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")
    add_resource_multiplier_arguments(database_parent)

    subparser_database = subparsers.add_parser("database", help="Install or update one managed annotation database")
    database_subparsers = subparser_database.add_subparsers(dest="database_name", help="Managed databases")
    database_subparsers.required = True

    database_kegg = database_subparsers.add_parser("kegg", parents=[database_parent], help="Install or update the KEGG/KOfam database", aliases=["kofams"])
    database_cazy = database_subparsers.add_parser("cazy", parents=[database_parent], help="Install or update the CAZy database")
    database_pfam = database_subparsers.add_parser("pfam", parents=[database_parent], help="Install or update the PFAM database")
    database_vfdb = database_subparsers.add_parser("vfdb", parents=[database_parent], help="Install or update the VFDB database")
    database_amr = database_subparsers.add_parser("amr", parents=[database_parent], help="Install or update the AMR database")

    subparser_environments = subparsers.add_parser("environments", help="Pre-create conda environments")
    subparser_environments.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_environments.add_argument("--profile", default="local", choices=["local", "slurm"])
    add_resource_multiplier_arguments(subparser_environments)

    subparser_unlock = subparsers.add_parser("unlock", help="Unlock snakemake")
    subparser_unlock.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_unlock.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")
    subparser_unlock.add_argument("-p", "--profile", required=False, default="slurm", help="Snakemake profile. Default is slurm")

    subparser_update = subparsers.add_parser("update", help="Reinstall Drakkar from the Git repo (forces reinstall in this environment)")
    subparser_update.add_argument("-e", "--env_path",type=str, help="Path to a shared conda environment directory (default: drakkar install path)")

    subparser_transfer = subparsers.add_parser("transfer", help="Transfer outputs via sftp")
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

    subparser_config = subparsers.add_parser("config", help="View or edit drakkar workflow/config.yaml")
    config_actions = subparser_config.add_mutually_exclusive_group(required=True)
    config_actions.add_argument("--view", action="store_true", help="Print workflow/config.yaml")
    config_actions.add_argument("--edit", action="store_true", help="Open workflow/config.yaml in a terminal editor")

    subparser_logging = subparsers.add_parser("logging", help="Inspect Drakkar run metadata and Snakemake logs")
    subparser_logging.add_argument("-o", "--output", required=False, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")
    subparser_logging.add_argument("--run", required=False, help="Specific run ID (YYYYMMDD-HHMMSS) or drakkar_<run_id>.yaml file name")
    subparser_logging.add_argument("--tail", required=False, type=positive_int, default=50, help="Number of log lines to show when no failure excerpt is found and --summary is not used. Default: 50")
    subparser_logging.add_argument("--summary", action="store_true", help="Print only the parsed workflow summary without log excerpts or tails")
    subparser_logging.add_argument("--full", action="store_true", help="Print the full Snakemake log instead of only the failure excerpt or tail")
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
            ("Start Here", ["complete"]),
            ("Data Generation Workflows", ["preprocessing", "cataloging"]),
            ("Analysis Workflows", ["profiling", "dereplicating", "annotating", "inspecting", "expressing"]),
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
            ("Workflow Scope", ["mode", "type", "annotation_type", "gtdb_version", "multicoverage", "fraction", "nonpareil", "ani"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Assembly Strategy", ["mode", "multicoverage"]),
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
            ("Run Configuration", ["output", "env_path", "profile", "overwrite"]),
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
                ("Run Configuration", ["env_path", "profile", "overwrite"]),
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
            ("Environment Setup", ["env_path", "profile"]),
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
        ],
        sections=[
            ("Execution", ["env_path"]),
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
            "drakkar logging -o drakkar_output --run 20260503-101530 --paths",
        ],
        sections=[
            ("Target Run", ["output", "run"]),
            ("Display Options", ["summary", "tail", "full", "paths", "list"]),
        ],
    )

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
        (getattr(args, "bins_file", None), "Bins file", False),
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
        return run_update()

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
            file_mags_to_json(bins_dir,args.output)
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
            file_mags_to_json(args.bins_dir,args.output)
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
            file_mags_to_json(args.bins_dir,args.output)
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


if __name__ == "__main__":
    main()
