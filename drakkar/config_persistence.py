"""Keep site-specific ``config.yaml`` values alive across Drakkar upgrades.

``workflow/config.yaml`` ships inside the package, so every reinstall — whether
it comes from ``drakkar update`` or a plain ``pip install --upgrade`` — replaces
it with the defaults baked into the new release. Everything a site configures
there is lost: the database releases installed by ``drakkar database``, the
paths ``--set-default`` wrote, and any hand-edited ``DATABASES_DIR`` or
``ENVIRONMENTS_DIR``.

Two things guard against that:

* every value worth keeping is mirrored into ``~/.drakkar/config-values.yaml``
  (``$DRAKKAR_HOME`` overrides the directory), which lives outside the package
  and therefore survives any reinstall, together with a timestamped copy of the
  whole file in ``~/.drakkar/config-backups/``;
* after an upgrade the saved values are reconciled against the freshly
  installed defaults and written back into ``config.yaml``.

Reconciliation never trusts a path blindly. For each key it looks at what
actually exists on disk:

* the saved path still exists — it is kept, unless a newer release of the same
  database sits next to it, in which case the newest one wins;
* the saved path is gone — the newest sibling release that does exist is used,
  falling back to the freshly installed default when that one exists;
* nothing exists — the saved value is kept anyway. A missing path on a login
  node usually means an unmounted filesystem, not a deleted database, so the
  configured value is never discarded on that evidence alone.

Keys that name one specific release on purpose (``GTDB_DB_226``) are only ever
preserved, never re-pointed at a newer release.
"""

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import yaml

from drakkar import __version__
from drakkar.cli_context import CONFIG_PATH
from drakkar.database_latest import _version_key as version_sort_key
from drakkar.database_registry import MANAGED_DATABASES
from drakkar.output import print

# Values worth carrying across a reinstall: the base directories and every
# database path. Module versions and workflow parameters are deliberately left
# out — those are what a new release is expected to change.
PRESERVED_KEY_RE = re.compile(r"^(?:DATABASES_DIR|ENVIRONMENTS_DIR|[A-Z0-9_]+_DB(?:_\d+)?|[A-Z0-9_]+_MODEL)$")

# GTDB_DB_226 and friends name one release on purpose: the annotating workflow
# resolves --gtdb-version to GTDB_DB_<version>, so re-pointing one at a newer
# release would silently answer a request for 226 with 232.
PINNED_KEY_RE = re.compile(r"^[A-Z0-9_]+_DB_\d+$")

# Keys whose value is a release directory sitting beside its other releases, so
# a newer one can be recognised. Bundled databases (foldseek) point at
# artifacts inside a release directory and are excluded.
RELEASE_KEYS = {
    definition["config_key"]
    for definition in MANAGED_DATABASES.values()
    if not definition.get("config_targets")
} | {"GTDB_DB"}

RELEASE_KEY_DATABASES = {
    definition["config_key"]: name
    for name, definition in MANAGED_DATABASES.items()
    if not definition.get("config_targets")
}

SOURCE_CONFIG = "kept"
SOURCE_INSTALLED = "installed default"
SOURCE_NEWEST = "newest on disk"

BACKUPS_KEPT = 10


@dataclass(frozen=True)
class Resolution:
    """What one config key ends up pointing at after an upgrade."""

    key: str
    value: str
    installed: str
    source: str
    detail: str = ""

    @property
    def changed(self) -> bool:
        """Whether the value differs from what the new release shipped."""
        return self.value != self.installed


def drakkar_home() -> Path:
    override = (os.environ.get("DRAKKAR_HOME") or "").strip()
    return Path(override).expanduser() if override else Path.home() / ".drakkar"


def store_path() -> Path:
    return drakkar_home() / "config-values.yaml"


def backup_dir() -> Path:
    return drakkar_home() / "config-backups"


def is_preserved_key(key) -> bool:
    return bool(PRESERVED_KEY_RE.match(str(key)))


def _clean(value) -> str:
    return "" if value is None else str(value).strip()


def read_config_values(config_path=None) -> dict:
    """The preserved keys currently written in config.yaml, empty ones included."""
    config_path = Path(config_path or CONFIG_PATH)
    if not config_path.is_file():
        return {}
    try:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return {}
    if not isinstance(config, dict):
        return {}
    return {key: _clean(value) for key, value in config.items() if is_preserved_key(key)}


def read_store() -> dict:
    """The values saved outside the package by an earlier install or update."""
    path = store_path()
    if not path.is_file():
        return {}
    try:
        saved = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return {}
    values = saved.get("values") if isinstance(saved, dict) else None
    if not isinstance(values, dict):
        return {}
    return {key: _clean(value) for key, value in values.items() if is_preserved_key(key) and _clean(value)}


def record_values(values) -> dict:
    """Merge non-empty values into the persistent store and return everything it holds."""
    merged = read_store()
    merged.update({key: _clean(value) for key, value in (values or {}).items() if is_preserved_key(key) and _clean(value)})
    path = store_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "saved_by": __version__,
        "saved_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "values": dict(sorted(merged.items())),
    }
    header = (
        "# Site-specific values Drakkar restores into workflow/config.yaml after an\n"
        "# upgrade replaces it. Written by 'drakkar database', 'drakkar update' and\n"
        "# 'drakkar config --restore'. Safe to edit or delete.\n"
    )
    path.write_text(header + yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return merged


def backup_config(config_path=None) -> Path | None:
    """Copy the whole config next to the store, keeping the last few copies."""
    config_path = Path(config_path or CONFIG_PATH)
    if not config_path.is_file():
        return None
    directory = backup_dir()
    directory.mkdir(parents=True, exist_ok=True)
    backup = directory / f"config-{datetime.now().strftime('%Y%m%d-%H%M%S')}.yaml"
    shutil.copy2(config_path, backup)
    stale = sorted(directory.glob("config-*.yaml"))[:-BACKUPS_KEPT]
    for path in stale:
        try:
            path.unlink()
        except OSError:
            pass
    return backup


def snapshot_config(config_path=None) -> dict:
    """Save the values in the current config before a reinstall replaces it.

    The config wins over the store wherever it holds something, so values that
    were hand-edited since the last save are the ones carried forward.
    """
    saved = read_store()
    saved.update({key: value for key, value in read_config_values(config_path).items() if value})
    backup_config(config_path)
    return record_values(saved)


def set_config_value(config_path, config_key, new_value) -> None:
    """Rewrite one key in place, leaving the rest of the file — comments included — untouched."""
    config_path = Path(config_path)
    pattern = re.compile(rf"^({re.escape(config_key)}:\s*)(\".*?\"|'.*?'|[^\n#]+)(\s*(#.*)?)?$", re.MULTILINE)
    config_text = config_path.read_text(encoding="utf-8")
    replacement = rf'\1"{new_value}"\3'
    updated_text, count = pattern.subn(replacement, config_text, count=1)
    if count != 1:
        raise ValueError(f"Could not update {config_key} in {config_path}")
    config_path.write_text(updated_text, encoding="utf-8")


def _exists(path) -> bool:
    try:
        return Path(path).exists()
    except OSError:
        return False


def _is_release_dir(key, path) -> bool:
    """Whether a directory beside the configured one looks like another release.

    Managed installs drop a ``database_versions.yaml`` next to the data, so an
    unrelated directory (a scratch or log directory) is not mistaken for one.
    """
    path = Path(path)
    if not path.is_dir() or not any(character.isdigit() for character in path.name):
        return False
    if (path / "database_versions.yaml").is_file():
        return True
    database = RELEASE_KEY_DATABASES.get(key)
    if database and (path / MANAGED_DATABASES[database]["basename"]).exists():
        return True
    # GTDB and other externally installed databases have no marker file, so the
    # release naming used by their own tooling is the only signal.
    return key not in RELEASE_KEY_DATABASES and bool(re.match(r"^release\d", path.name))


def _sibling_releases(key, configured) -> list[str]:
    parent = Path(configured).parent
    if not parent.is_dir():
        return []
    try:
        entries = sorted(parent.iterdir())
    except OSError:
        return []
    return [str(entry) for entry in entries if _is_release_dir(key, entry)]


def resolve_value(key, saved, installed):
    """Pick the path one key should hold: the saved one, the newest, or the shipped default."""
    saved = _clean(saved)
    installed = _clean(installed)
    if not saved:
        return installed, SOURCE_INSTALLED, ""
    if saved == installed:
        return saved, SOURCE_CONFIG, ""
    if PINNED_KEY_RE.match(key):
        detail = "" if _exists(saved) else "path not found; kept as configured"
        return saved, SOURCE_CONFIG, detail

    candidates = {}
    if _exists(saved):
        candidates[saved] = SOURCE_CONFIG
    if installed and _exists(installed):
        candidates.setdefault(installed, SOURCE_INSTALLED)
    if key in RELEASE_KEYS:
        for release in _sibling_releases(key, saved) + (_sibling_releases(key, installed) if installed else []):
            candidates.setdefault(release, SOURCE_NEWEST)
    if not candidates:
        return saved, SOURCE_CONFIG, "no release found on disk; kept as configured"

    # Newest release wins; a tie keeps the configured path.
    newest = max(candidates, key=lambda path: (version_sort_key(Path(path).name), path == saved))
    return newest, candidates[newest], ""


def reconcile_config(config_path=None, saved=None) -> list[Resolution]:
    """Write the saved values back into a freshly installed config.yaml."""
    config_path = Path(config_path or CONFIG_PATH)
    saved = read_store() if saved is None else saved
    installed = read_config_values(config_path)
    resolutions = []
    for key, installed_value in installed.items():
        value, source, detail = resolve_value(key, saved.get(key), installed_value)
        resolution = Resolution(key=key, value=value, installed=installed_value, source=source, detail=detail)
        if resolution.changed:
            try:
                set_config_value(config_path, key, value)
            except (OSError, ValueError) as error:
                resolution = Resolution(
                    key=key,
                    value=installed_value,
                    installed=installed_value,
                    source=SOURCE_INSTALLED,
                    detail=f"could not be restored: {error}",
                )
        resolutions.append(resolution)
    record_values({resolution.key: resolution.value for resolution in resolutions})
    return resolutions


def print_reconcile_report(resolutions, config_path=None) -> None:
    """Summarise what the reconciliation kept, upgraded, and could not find."""
    config_path = Path(config_path or CONFIG_PATH)
    if not resolutions:
        print(f"No database values were found in {config_path}.")
        return

    restored = [item for item in resolutions if item.changed]
    failed = [item for item in resolutions if item.detail and not item.changed]

    if restored:
        print(f"Restored {len(restored)} database value(s) into {config_path}:")
        for item in restored:
            print(f"  {item.key}: {item.value}")
            if item.source == SOURCE_NEWEST:
                print("    re-pointed at the newest release found on disk")
            elif item.detail:
                print(f"    {item.detail}")
    else:
        print(f"Database configuration in {config_path.name} already matches the saved values.")

    for item in failed:
        print(f"  {item.key}: {item.detail}")

    print(f"Saved values are kept in {store_path()}.")
