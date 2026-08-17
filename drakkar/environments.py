"""Inspection and cleanup of the shared conda environment directory.

Snakemake deploys every ``conda:`` environment into ``--conda-prefix`` as a
directory named after a hash of the environment definition, next to a verbatim
copy of that definition (``<hash>.yaml``) and, when present, its post-deploy
script (``<hash>.post-deploy.sh``). Because the hash covers the file content,
editing any ``workflow/envs/*.yaml`` produces a brand new directory and leaves
the previous one behind forever.

Rather than reimplementing Snakemake's hashing (which changes between
releases), environments are matched by comparing the sidecar copy against the
environment definitions shipped by the installed Drakkar version.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from datetime import datetime
from pathlib import Path

import yaml

from drakkar import __version__
from drakkar.cli_context import ERROR, INFO, PACKAGE_DIR, RESET
from drakkar.output import print

ENVS_DIR = PACKAGE_DIR / "workflow" / "envs"

# Snakemake names deployed environments after a hex digest. Anything else in
# the directory is left alone: a shared prefix may hold unrelated environments.
HASH_NAME_RE = re.compile(r"^[0-9a-f]{6,64}$")

ENV_SUFFIXES = (".yaml", ".yml")
SIDECAR_SUFFIXES = (".yaml", ".yml", ".post-deploy.sh", ".pin.txt")

STATUS_IN_USE = "in use"
STATUS_ORPHAN = "orphan"
STATUS_INCOMPLETE = "incomplete"
STATUS_UNKNOWN = "unknown"

REMOVABLE_STATUSES = (STATUS_ORPHAN, STATUS_INCOMPLETE)

# Report order: what is live first, then what can be reclaimed, then what is
# left untouched.
STATUS_ORDER = (STATUS_IN_USE, STATUS_ORPHAN, STATUS_INCOMPLETE, STATUS_UNKNOWN)

STATUS_LABELS = {
    STATUS_ORPHAN: "(superseded definition)",
    STATUS_INCOMPLETE: "(no definition file)",
    STATUS_UNKNOWN: "(not a Drakkar environment)",
}


def fingerprint(text):
    """Return a comparable form of an environment definition.

    The YAML structure is compared rather than raw bytes so that whitespace or
    key ordering differences do not mark a live environment as an orphan. List
    order is preserved because channel order is meaningful.
    """
    try:
        parsed = yaml.safe_load(text)
    except yaml.YAMLError:
        parsed = None
    if parsed is None:
        return "\n".join(line.rstrip() for line in text.splitlines() if line.strip())
    return json.dumps(parsed, sort_keys=True, default=str)


def package_environments(envs_dir=None):
    """Map fingerprint -> environment file name for the installed version."""
    envs_path = Path(envs_dir) if envs_dir else ENVS_DIR
    known = {}
    if not envs_path.is_dir():
        return known
    for env_file in sorted(envs_path.iterdir()):
        if not env_file.is_file() or env_file.suffix not in ENV_SUFFIXES:
            continue
        try:
            content = env_file.read_text(encoding="utf-8")
        except OSError:
            continue
        known.setdefault(fingerprint(content), env_file.name)
    return known


def directory_size(path):
    """Total size of a directory tree in bytes, ignoring unreadable entries."""
    total = 0
    for root, dirs, files in os.walk(path, onerror=lambda error: None):
        for name in files:
            try:
                stat = os.lstat(os.path.join(root, name))
            except OSError:
                continue
            total += stat.st_size
    return total


def format_size(size):
    if size is None:
        return "-"
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            precision = 0 if unit == "B" else 1
            return f"{value:.{precision}f} {unit}"
        value /= 1024
    return f"{value:.1f} TiB"


def format_date(timestamp):
    if timestamp is None:
        return "-"
    return datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d")


def _sidecar_paths(env_dir, env_hash):
    return [
        candidate
        for suffix in SIDECAR_SUFFIXES
        for candidate in [env_dir / f"{env_hash}{suffix}"]
        if candidate.exists()
    ]


def _sidecar_definition(sidecars):
    for path in sidecars:
        if path.suffix in ENV_SUFFIXES:
            try:
                return path.read_text(encoding="utf-8")
            except OSError:
                return None
    return None


def scan_environments(env_path, compute_size=True, envs_dir=None):
    """Classify every entry of the conda prefix directory.

    Returns a dict with the resolved directory, the classified entries, the
    environment definitions of the installed version that are not deployed yet,
    and any leftover sidecar files whose environment directory is gone.
    """
    env_dir = Path(env_path).expanduser()
    known = package_environments(envs_dir)

    entries = []
    seen_hashes = set()
    deployed_names = set()

    for child in sorted(env_dir.iterdir(), key=lambda item: item.name):
        if not child.is_dir():
            continue
        env_hash = child.name
        seen_hashes.add(env_hash)

        if not HASH_NAME_RE.match(env_hash):
            entries.append(
                {
                    "hash": env_hash,
                    "path": child,
                    "name": None,
                    "status": STATUS_UNKNOWN,
                    "sidecars": [],
                    "size": directory_size(child) if compute_size else None,
                    "created": _mtime(child),
                    "removable": False,
                }
            )
            continue

        sidecars = _sidecar_paths(env_dir, env_hash)
        definition = _sidecar_definition(sidecars)
        if definition is None:
            status = STATUS_INCOMPLETE
            name = None
        else:
            name = known.get(fingerprint(definition))
            status = STATUS_IN_USE if name else STATUS_ORPHAN
            if name:
                deployed_names.add(name)

        entries.append(
            {
                "hash": env_hash,
                "path": child,
                "name": name,
                "status": status,
                "sidecars": sidecars,
                "size": directory_size(child) if compute_size else None,
                "created": _mtime(child),
                "removable": status in REMOVABLE_STATUSES,
            }
        )

    leftovers = []
    for child in sorted(env_dir.iterdir(), key=lambda item: item.name):
        if child.is_dir():
            continue
        env_hash = _sidecar_hash(child.name)
        if env_hash is None or env_hash in seen_hashes:
            continue
        if not HASH_NAME_RE.match(env_hash):
            continue
        leftovers.append(child)

    missing = sorted(name for name in set(known.values()) if name not in deployed_names)
    entries.sort(key=lambda entry: (STATUS_ORDER.index(entry["status"]), entry["name"] or "", entry["hash"]))

    return {
        "directory": env_dir,
        "entries": entries,
        "leftovers": leftovers,
        "missing": missing,
        "known": known,
    }


def _mtime(path):
    try:
        return path.stat().st_mtime
    except OSError:
        return None


def _sidecar_hash(filename):
    for suffix in SIDECAR_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return None


def _validate_directory(env_path):
    if not env_path:
        print(f"{ERROR}ERROR:{RESET} No conda environment directory configured.")
        return None
    env_dir = Path(env_path).expanduser()
    if not env_dir.exists():
        print(f"{ERROR}ERROR:{RESET} Environment directory not found: {env_dir}")
        return None
    if not env_dir.is_dir():
        print(f"{ERROR}ERROR:{RESET} Environment path is not a directory: {env_dir}")
        return None
    return env_dir


def _print_table(entries, compute_size):
    if not entries:
        print("No conda environments found in this directory.")
        return
    print(f"{'ENVIRONMENT':<32} {'HASH':<14} {'SIZE':>9}  {'CREATED':<12} STATUS")
    for entry in entries:
        name = entry["name"] or STATUS_LABELS.get(entry["status"], "(unknown)")
        size = format_size(entry["size"]) if compute_size else "-"
        print(
            f"{name[:32]:<32} {entry['hash'][:14]:<14} {size:>9}  "
            f"{format_date(entry['created']):<12} {entry['status']}"
        )


def _summarise(scan, compute_size, include_missing=True):
    entries = scan["entries"]
    removable = [entry for entry in entries if entry["removable"]]
    in_use = [entry for entry in entries if entry["status"] == STATUS_IN_USE]
    unknown = [entry for entry in entries if entry["status"] == STATUS_UNKNOWN]

    print("")
    print(
        f"{len(entries)} environment(s): {len(in_use)} in use, "
        f"{len(removable)} reclaimable, {len(unknown)} unknown."
    )
    if removable and compute_size:
        reclaimable = sum(entry["size"] or 0 for entry in removable)
        print(f"Reclaimable space: {format_size(reclaimable)}")
    if scan["leftovers"]:
        print(f"Leftover definition files without an environment: {len(scan['leftovers'])}")
    if include_missing and scan["missing"]:
        print(f"{INFO}INFO:{RESET} Not built yet: {', '.join(scan['missing'])}")
    return removable


def run_environments_list(env_path, compute_size=True):
    env_dir = _validate_directory(env_path)
    if env_dir is None:
        return 1

    scan = scan_environments(env_dir, compute_size=compute_size)
    print(f"Environment directory: {scan['directory']}")
    print(f"Drakkar definitions: {len(set(scan['known'].values()))} in {ENVS_DIR}")
    print("")
    _print_table(scan["entries"], compute_size)
    _summarise(scan, compute_size)

    unknown = [entry for entry in scan["entries"] if entry["status"] == STATUS_UNKNOWN]
    if unknown:
        print("")
        print(
            f"{INFO}INFO:{RESET} Directories not named like a Snakemake environment "
            "are reported as unknown and are never removed by --prune."
        )
    return 0


def run_environments_prune(env_path, assume_yes=False, compute_size=True):
    env_dir = _validate_directory(env_path)
    if env_dir is None:
        return 1

    scan = scan_environments(env_dir, compute_size=compute_size)
    removable = [entry for entry in scan["entries"] if entry["removable"]]
    leftovers = scan["leftovers"]

    print(f"Environment directory: {scan['directory']}")
    print("")
    if not removable and not leftovers:
        print("Nothing to prune: every environment matches the installed Drakkar version.")
        return 0

    _print_table(removable, compute_size)
    for leftover in leftovers:
        print(f"{'(leftover file)':<32} {leftover.name}")
    _summarise(scan, compute_size, include_missing=False)

    print("")
    print(
        f"{INFO}INFO:{RESET} Environments are matched against Drakkar {__version__}. "
        "If another Drakkar version shares this directory, the environments it "
        "needs are reported here as reclaimable."
    )

    if not assume_yes:
        print("")
        print(
            f"{INFO}INFO:{RESET} Dry run. Rerun with --prune --yes to delete the "
            "environments listed above."
        )
        return 0

    print("")
    removed = 0
    failed = 0
    for entry in removable:
        targets = [entry["path"], *entry["sidecars"]]
        if _remove_paths(targets, entry["hash"]):
            removed += 1
        else:
            failed += 1
    for leftover in leftovers:
        if _remove_paths([leftover], leftover.name):
            removed += 1
        else:
            failed += 1

    print("")
    print(f"Removed {removed} item(s) from {scan['directory']}.")
    if failed:
        print(f"{ERROR}ERROR:{RESET} {failed} item(s) could not be removed.")
        return 1
    return 0


def _remove_paths(paths, label):
    for path in paths:
        try:
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()
        except OSError as error:
            print(f"{ERROR}ERROR:{RESET} Could not remove {path}: {error}")
            return False
    print(f"Removed {label}")
    return True
