import os
import re
import shlex
import shutil
import subprocess
import sys

from drakkar.cli_context import CONFIG_PATH, ERROR, RESET
from drakkar.database_registry import MANAGED_DATABASES, database_release_dir
from drakkar.output import print

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
