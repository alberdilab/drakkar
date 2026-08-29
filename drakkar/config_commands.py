import os
import shlex
import shutil
import subprocess
import sys

from drakkar.cli_context import CONFIG_PATH, ERROR, INFO, RESET
from drakkar.config_persistence import (
    backup_config,
    print_reconcile_report,
    read_store,
    reconcile_config,
    record_values,
    set_config_value,
    store_path,
)
from drakkar.database_registry import database_config_targets
from drakkar.output import print

def replace_config_value(config_key, new_value):
    set_config_value(CONFIG_PATH, config_key, new_value)

def set_default_database_path(database_name, directory, version):
    # Returns {config_key: path} for every key this database sets. Single-target
    # databases yield one entry; bundled ones (e.g. foldseek) yield several.
    targets = database_config_targets(database_name, directory, version)
    for config_key, path in targets.items():
        replace_config_value(config_key, path)
    # Mirrored outside the package so the next reinstall, which replaces
    # config.yaml wholesale, can put these paths back.
    record_values(targets)
    return targets

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

def restore_config():
    """Put the saved database values back into a config.yaml an upgrade replaced."""
    if not CONFIG_PATH.exists():
        print(f"{ERROR}ERROR:{RESET} config.yaml not found: {CONFIG_PATH}")
        return 1
    if not read_store():
        print(f"{ERROR}ERROR:{RESET} No saved configuration values found: {store_path()}")
        print(f"{INFO}Values are saved by 'drakkar database', 'drakkar update' and this command.{RESET}")
        return 1
    backup_config(CONFIG_PATH)
    print_reconcile_report(reconcile_config(CONFIG_PATH), CONFIG_PATH)
    return 0
