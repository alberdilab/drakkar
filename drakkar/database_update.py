"""Install the newest release of every managed database that has fallen behind.

``drakkar database update`` is ``drakkar database latest`` followed by the
install command it prints, for each database at once. It reuses the discovery in
:mod:`drakkar.database_latest` to decide what is outdated, then runs the ordinary
per-database install workflow once per database.

The installs run one after another rather than as a single Snakemake session on
purpose. Each release directory *is* the Snakemake working directory for its
install, which is what makes ``drakkar logging -o <release_dir>``,
``drakkar unlock -o <release_dir>``, the per-release failure report and
``database_versions.yaml`` work. Running them sequentially keeps every one of
those semantics untouched and isolates failures: one database that cannot be
downloaded leaves the others installed and reports where its own failure is
recorded.

Because an update downloads tens of gigabytes and then repoints ``config.yaml``
at the new releases, the command prints the plan and stops unless ``--yes`` is
given.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from drakkar.cli_context import ERROR, INFO, RESET
from drakkar.database_latest import (
    DEFAULT_TIMEOUT,
    LATEST_SOURCES,
    STATUS_AHEAD,
    STATUS_CURRENT,
    STATUS_OUTDATED,
    normalize_latest_names,
    resolve_all,
)
from drakkar.database_registry import (
    MANAGED_DATABASES,
    database_release_dir,
    database_release_from_config,
)
from drakkar.output import print


@dataclass(frozen=True)
class DatabaseUpdate:
    """One database to install, and where its new release will land."""

    name: str
    label: str
    config_key: str
    installed: str
    latest: str
    base_directory: Path
    release_dir: Path


def _base_directory(name, configured):
    """The directory the release folders of this database live in."""
    return database_release_from_config(name, Path(configured)).parent


def plan_database_updates(names=None, timeout=DEFAULT_TIMEOUT, config=None):
    """Split the checked databases into the ones to install and the ones to skip.

    Returns ``(updates, skipped)``, where ``skipped`` pairs each untouched
    database with the reason it is being left alone.
    """
    results = resolve_all(names, config=config, timeout=timeout)
    updates, skipped = [], []
    for result in results:
        if result.status == STATUS_CURRENT:
            continue
        if result.status == STATUS_AHEAD:
            skipped.append((result, f"configured release {result.installed} is newer than the source lists"))
            continue
        if result.status != STATUS_OUTDATED:
            skipped.append((result, result.detail or "could not be checked"))
            continue
        if result.name not in MANAGED_DATABASES:
            hint = LATEST_SOURCES[result.name].get("install_hint") or "not installed by Drakkar"
            skipped.append((result, f"{result.latest} is available, but this database is {hint}"))
            continue
        base_directory = _base_directory(result.name, result.configured)
        updates.append(
            DatabaseUpdate(
                name=result.name,
                label=result.label,
                config_key=result.config_key,
                installed=result.installed,
                latest=result.latest,
                base_directory=base_directory,
                release_dir=database_release_dir(result.name, base_directory, result.latest),
            )
        )
    return updates, skipped


def print_update_plan(updates, skipped, checked, assume_yes):
    """Show what would be installed before anything is downloaded."""
    if updates:
        print("The following databases will be installed:")
        print("")
        for update in updates:
            print(f"  {update.name} {update.installed} -> {update.latest}")
            print(f"    into: {update.release_dir}")
            print(f"    config.yaml key: {update.config_key}")
        print("")
    for result, reason in skipped:
        print(f"{INFO}INFO:{RESET} Skipping {result.name}: {reason}")
    if skipped:
        print("")

    up_to_date = checked - len(updates) - len(skipped)
    print(f"{checked} database(s) checked: {len(updates)} to install, {up_to_date} up to date, {len(skipped)} skipped.")

    if not updates:
        return
    if not assume_yes:
        print("")
        print(f"{INFO}INFO:{RESET} Dry run. Rerun with --yes to download and install the releases listed above.")


def apply_database_updates(updates, install, set_default=True):
    """Install each release in turn, keeping going when one of them fails."""
    installed, failed = [], []
    for update in updates:
        print("")
        print(f"{INFO}INFO:{RESET} Installing {update.name} {update.latest} into {update.release_dir}")
        if install(update):
            installed.append(update)
        else:
            failed.append(update)
            print(f"{ERROR}ERROR:{RESET} {update.name} was not installed. The releases already installed are kept.")

    print("")
    print(f"{len(updates)} database(s): {len(installed)} installed, {len(failed)} failed.")
    if failed:
        print("")
        print(f"{ERROR}ERROR:{RESET} These databases were not updated: {', '.join(update.name for update in failed)}")
        for update in failed:
            print(f"  inspect with: drakkar logging -o {update.release_dir} --failures")
    if installed and not set_default:
        print("")
        print(f"{INFO}INFO:{RESET} config.yaml was left unchanged (--no-set-default). Point these keys at the new releases to use them:")
        for update in installed:
            print(f"  {update.config_key}: {update.release_dir}")
    if installed and set_default:
        print("")
        print(
            f"{INFO}INFO:{RESET} Existing outputs were built with the earlier releases. "
            "Rerunning a directory that holds them needs --allow-database-change, or "
            "the affected outputs deleted so they are rebuilt."
        )
    return 1 if failed else 0


def run_database_update(
    names=None,
    install=None,
    timeout=DEFAULT_TIMEOUT,
    assume_yes=False,
    set_default=True,
    config=None,
):
    """Check every database and install the releases that are behind."""
    selected, unknown = normalize_latest_names(names)
    if unknown:
        print(f"{ERROR}ERROR:{RESET} Unknown database(s): {', '.join(str(value) for value in unknown)}")
        print(f"Databases that can be checked: {', '.join(LATEST_SOURCES)}")
        return 1

    print("Checking the configured databases against their sources.")
    print("")
    updates, skipped = plan_database_updates(selected, timeout=timeout, config=config)
    print_update_plan(updates, skipped, len(selected), assume_yes)

    if not updates or not assume_yes:
        return 0
    return apply_database_updates(updates, install, set_default=set_default)
