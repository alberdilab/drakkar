"""Report the newest release each database publishes upstream.

``drakkar database latest`` answers one question: is the release wired into
``config.yaml`` still the newest one the source offers? Databases are installed
once and then quietly go stale for months, and nothing in a workflow run
surfaces that — the provenance checks in :mod:`drakkar.database_checks` only
compare a run against earlier runs in the same output directory, never against
the upstream source.

Sources expose their releases in three different ways, so each database
declares which one applies:

* ``index`` — the source publishes a browsable page whose links identify the
  releases (KEGG/KOfam, Pfam, NCBIfam-AMRFinder, AMRFinderPlus, CARD).
* ``probe`` — the source has no listing at all, so consecutive version numbers
  are requested until one is missing (dbCAN/CAZy).
* ``mtime`` — the source is a single rolling file with no versions, so its
  ``Last-Modified`` date is compared against the download date the release
  directory was named after (VFDB, the UniProt half of Foldseek).
* ``text`` — the source publishes its current version in a small text file
  (GTDB).

Every probe is a network call to a third-party host, so a failure is reported
as ``unknown`` for that one database and never fails the command: a stale
mirror must not stop someone from checking the other six.
"""

from __future__ import annotations

import re
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

import yaml

from drakkar import __version__
from drakkar.cli_context import ERROR, INFO, RESET, config_vars
from drakkar.database_registry import (
    MANAGED_DATABASES,
    database_base_directory,
    database_release_from_config,
    normalize_managed_database_name,
)
from drakkar.output import print

DEFAULT_TIMEOUT = 20
USER_AGENT = f"drakkar/{__version__}"

STATUS_CURRENT = "up to date"
STATUS_OUTDATED = "outdated"
STATUS_AHEAD = "ahead of source"
STATUS_UNKNOWN = "unknown"
STATUS_UNCONFIGURED = "not configured"

STATUS_ORDER = (
    STATUS_OUTDATED,
    STATUS_AHEAD,
    STATUS_UNKNOWN,
    STATUS_UNCONFIGURED,
    STATUS_CURRENT,
)

# How to discover the newest release of each database, and how to read the
# installed version out of the path stored in config.yaml.
#
# ``managed`` databases are the ones ``drakkar database <name>`` installs: their
# release directory is named after the version, so the installed version is the
# directory name. The others are installed by their own tooling, so a pattern
# pulls the version out of whatever the configured path happens to look like.
LATEST_SOURCES = {
    "kegg": {
        "label": "KEGG/KOfam",
        "config_key": "KEGG_DB",
        "managed": True,
        "strategy": "index",
        "url": "https://www.genome.jp/ftp/db/kofam/archives/",
        "entry_pattern": r"^(\d{4}-\d{2}-\d{2})/$",
    },
    "cazy": {
        "label": "CAZy/dbCAN",
        "config_key": "CAZY_DB",
        "managed": True,
        # dbCAN serves releases through a PHP endpoint with no directory
        # listing, so the next version is probed directly. The endpoint ignores
        # Range and streams the whole ~120 MB file on GET, so this must stay a
        # HEAD request.
        "strategy": "probe",
        "url_template": "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/{version}/dbCAN-HMMdb-{version}.txt",
        "version_template": "V{number}",
        "installed_pattern": r"^(V\d+)$",
        "probe_floor": 14,
        "probe_limit": 20,
    },
    "pfam": {
        "label": "Pfam",
        "config_key": "PFAM_DB",
        "managed": True,
        # Pfam_SARS-CoV-2_* releases share the directory and are excluded by
        # requiring the plain PfamNN.N form.
        "strategy": "index",
        "url": "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/",
        "entry_pattern": r"^(Pfam\d+\.\d+)/$",
    },
    "vfdb": {
        "label": "VFDB",
        "config_key": "VFDB_DB",
        "managed": True,
        "strategy": "mtime",
        "url": "https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz",
        "installed_pattern": r"^(\d{4}-\d{2}-\d{2})$",
    },
    "amr": {
        "label": "NCBIfam-AMRFinder",
        "config_key": "AMR_DB",
        "managed": True,
        "strategy": "index",
        "url": "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/",
        "entry_pattern": r"^(\d{4}-\d{2}-\d{2}\.\d+)/$",
    },
    "amrfinderplus": {
        "label": "AMRFinderPlus database",
        "config_key": "AMRFINDER_DB",
        "managed": True,
        "strategy": "index",
        # The AMR workflow pins AMRFinderPlus 4.2.x. NCBI partitions database
        # releases by compatible software minor version, so this deliberately
        # reports the newest database the bundled executable can consume.
        "url": "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/",
        "entry_pattern": r"^(\d{4}-\d{2}-\d{2}\.\d+)/$",
    },
    "card": {
        "label": "CARD",
        "config_key": "CARD_DB",
        "managed": True,
        "strategy": "index",
        "url": "https://card.mcmaster.ca/download",
        "entry_pattern": r"^/download/0/broadstreet-v(\d+\.\d+\.\d+)\.tar\.bz2$",
    },
    "foldseek": {
        "label": "Foldseek/UniProt Swiss-Prot",
        "config_key": "FOLDSEEK_DB",
        "managed": True,
        # The AlphaFold and ProstT5 halves are unversioned downloads; UniProt
        # Swiss-Prot is the only dated component of the bundle.
        "strategy": "mtime",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
        "installed_pattern": r"^(\d{4}-\d{2}-\d{2})$",
    },
    "gtdb": {
        "label": "GTDB-Tk reference data",
        "config_key": "GTDB_DB",
        "managed": False,
        "strategy": "text",
        "url": "https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt",
        "value_pattern": r"v(\d+(?:\.\d+)?)",
        "installed_pattern": r"^release(\d+(?:\.\d+)?)$",
        "install_hint": "install with GTDB-Tk (gtdbtk download-db.sh) and update GTDB_DB in config.yaml",
    },
}


@dataclass(frozen=True)
class LatestRelease:
    """The installed and newest-available release of one database."""

    name: str
    label: str
    config_key: str
    configured: str | None
    installed: str | None
    latest: str | None
    status: str
    detail: str = ""
    install_command: str | None = None


def _version_key(value):
    """Natural sort key: Pfam38.2 sorts above Pfam37.4, 2026-07-02 above 2026-02-01.

    Every element is a same-shaped triple so keys built from differently shaped
    version strings still compare without raising.
    """
    parts = [part for part in re.split(r"(\d+)", str(value)) if part]
    return tuple(
        (0, int(part), "") if part.isdigit() else (1, 0, part.lower())
        for part in parts
    )


def _open(url, method="GET", timeout=DEFAULT_TIMEOUT):
    request = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def _fetch_text(url, timeout=DEFAULT_TIMEOUT, limit=4_000_000):
    with _open(url, timeout=timeout) as response:
        return response.read(limit).decode("utf-8", errors="replace")


def _latest_from_index(definition, timeout):
    """Newest release directory listed in a browsable source directory."""
    html = _fetch_text(definition["url"], timeout=timeout)
    pattern = re.compile(definition["entry_pattern"])
    versions = {
        match.group(1)
        for href in re.findall(r'href="([^"]+)"', html)
        if (match := pattern.match(href))
    }
    if not versions:
        raise ValueError(f"no releases matched in {definition['url']}")
    return max(versions, key=_version_key), ""


def _url_exists(url, timeout):
    try:
        with _open(url, method="HEAD", timeout=timeout):
            return True
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        raise


def _latest_from_probe(definition, installed, timeout):
    """Newest release of a source with no listing, found by probing forward."""
    match = re.search(r"(\d+)", installed or "")
    number = int(match.group(1)) if match else definition["probe_floor"]
    if not _url_exists(definition["url_template"].format(version=definition["version_template"].format(number=number)), timeout):
        # The installed version is not offered any more; fall back to the floor
        # so the probe still reports something meaningful.
        number = definition["probe_floor"]
    newest = number
    for candidate in range(number + 1, number + definition["probe_limit"] + 1):
        version = definition["version_template"].format(number=candidate)
        if not _url_exists(definition["url_template"].format(version=version), timeout):
            break
        newest = candidate
    return definition["version_template"].format(number=newest), ""


def _latest_from_mtime(definition, timeout):
    """Publication date of a rolling, unversioned download."""
    with _open(definition["url"], method="HEAD", timeout=timeout) as response:
        last_modified = response.headers.get("Last-Modified")
    if not last_modified:
        raise ValueError(f"no Last-Modified header on {definition['url']}")
    published = parsedate_to_datetime(last_modified).astimezone(timezone.utc)
    return published.strftime("%Y-%m-%d"), "rolling release, compared by source date"


def _latest_from_text(definition, timeout):
    """Version published in a small text file by the source."""
    text = _fetch_text(definition["url"], timeout=timeout, limit=4096)
    match = re.search(definition["value_pattern"], text)
    if not match:
        raise ValueError(f"no version found in {definition['url']}")
    return match.group(1), ""


def configured_value(config, config_key):
    value = (config or {}).get(config_key)
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def installed_version(name, definition, configured):
    """The release currently wired into config.yaml, as a comparable version.

    Managed installs record the version they were asked for in
    ``database_versions.yaml``; that record wins over the directory name, which
    a user is free to rename.
    """
    if not configured:
        return None
    path = Path(configured)
    if definition.get("managed") and name in MANAGED_DATABASES:
        path = database_release_from_config(name, path)
        recorded = _recorded_version(path)
        if recorded:
            return recorded
    version = path.name
    pattern = definition.get("installed_pattern")
    if not pattern:
        return version or None
    match = re.match(pattern, version)
    return match.group(1) if match else None


def _recorded_version(release_dir):
    versions_file = Path(release_dir) / "database_versions.yaml"
    if not versions_file.is_file():
        return None
    try:
        recorded = yaml.safe_load(versions_file.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return None
    value = recorded.get("requested_version")
    return str(value) if value else None


def _install_command(name, definition, configured, latest, config=None):
    if not definition.get("managed") or not latest:
        return definition.get("install_hint")
    if configured:
        directory = str(database_release_from_config(name, Path(configured)).parent)
    else:
        # Nothing configured to install beside, so fall back to the base
        # directory the install commands default to.
        base = database_base_directory(name, (config or {}).get("DATABASES_DIR"))
        directory = str(base) if base else "<directory>"
    # Rolling sources have no version to ask for: the installer names the
    # release directory after the date it downloaded them.
    if definition["strategy"] == "mtime":
        return f"drakkar database {name} --directory {directory}"
    return f"drakkar database {name} --directory {directory} --version {latest}"


def _compare(installed, latest):
    if installed == latest:
        return STATUS_CURRENT
    return STATUS_OUTDATED if _version_key(installed) < _version_key(latest) else STATUS_AHEAD


def resolve_latest(name, config=None, timeout=DEFAULT_TIMEOUT):
    """Look up one database upstream and compare it with the configured release."""
    definition = LATEST_SOURCES[name]
    config = config_vars if config is None else config
    configured = configured_value(config, definition["config_key"])
    installed = installed_version(name, definition, configured)

    def result(latest, status, detail=""):
        return LatestRelease(
            name=name,
            label=definition["label"],
            config_key=definition["config_key"],
            configured=configured,
            installed=installed,
            latest=latest,
            status=status,
            detail=detail,
            install_command=_install_command(name, definition, configured, latest, config),
        )

    try:
        strategy = definition["strategy"]
        if strategy == "index":
            latest, detail = _latest_from_index(definition, timeout)
        elif strategy == "probe":
            latest, detail = _latest_from_probe(definition, installed, timeout)
        elif strategy == "mtime":
            latest, detail = _latest_from_mtime(definition, timeout)
        else:
            latest, detail = _latest_from_text(definition, timeout)
    except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError) as error:
        return result(None, STATUS_UNKNOWN, f"could not reach the source: {error}")

    if not configured:
        return result(latest, STATUS_UNCONFIGURED, f"{definition['config_key']} is empty in config.yaml")
    if not installed:
        return result(latest, STATUS_UNKNOWN, f"no release version could be read from the configured path: {configured}")
    return result(latest, _compare(installed, latest), detail)


def resolve_all(names=None, config=None, timeout=DEFAULT_TIMEOUT):
    """Query every requested database, one thread each so probes overlap."""
    selected = list(names or LATEST_SOURCES)
    if not selected:
        return []
    with ThreadPoolExecutor(max_workers=len(selected)) as pool:
        results = list(pool.map(lambda name: resolve_latest(name, config=config, timeout=timeout), selected))
    return sorted(results, key=lambda item: (STATUS_ORDER.index(item.status), item.name))


def normalize_latest_names(values):
    """Resolve user-supplied names (and managed aliases) to source names."""
    if not values:
        return list(LATEST_SOURCES), []
    selected, unknown = [], []
    for value in values:
        name = str(value or "").strip().lower()
        if name not in LATEST_SOURCES:
            # kofams -> kegg, so the check accepts the same names as the installer.
            name = normalize_managed_database_name(name) or name
        if name in LATEST_SOURCES:
            if name not in selected:
                selected.append(name)
        else:
            unknown.append(value)
    return selected, unknown


def _print_table(results):
    print(f"{'DATABASE':<22} {'CONFIGURED':<16} {'LATEST':<16} STATUS")
    for result in results:
        installed = result.installed or "-"
        latest = result.latest or "-"
        print(f"{result.name[:22]:<22} {installed[:16]:<16} {latest[:16]:<16} {result.status}")


def run_database_latest(names=None, timeout=DEFAULT_TIMEOUT, config=None):
    """Report, for each database, whether config.yaml still points at the newest release."""
    selected, unknown = normalize_latest_names(names)
    if unknown:
        print(f"{ERROR}ERROR:{RESET} Unknown database(s): {', '.join(str(value) for value in unknown)}")
        print(f"Databases that can be checked: {', '.join(LATEST_SOURCES)}")
        return 1

    print(f"Checked on {datetime.now(timezone.utc).strftime('%Y-%m-%d')} (UTC) against the upstream sources.")
    print("")
    results = resolve_all(selected, config=config, timeout=timeout)
    _print_table(results)

    def with_status(*statuses):
        return [result for result in results if result.status in statuses]

    outdated = with_status(STATUS_OUTDATED)
    ahead = with_status(STATUS_AHEAD)
    unresolved = with_status(STATUS_UNKNOWN, STATUS_UNCONFIGURED)

    for result in outdated:
        print("")
        print(f"{result.label} ({result.config_key}): {result.installed} -> {result.latest}")
        if result.detail:
            print(f"  {result.detail}")
        if result.install_command:
            print(f"  update with: {result.install_command}")

    for result in ahead:
        print("")
        print(
            f"{result.label} ({result.config_key}): {result.installed} is newer than "
            f"the newest release the source lists ({result.latest})."
        )

    for result in unresolved:
        print("")
        print(f"{result.label} ({result.config_key}): {result.detail}")

    print("")
    print(
        f"{len(results)} database(s): {len(with_status(STATUS_CURRENT))} up to date, "
        f"{len(outdated)} outdated, {len(ahead)} ahead of source, {len(unresolved)} unresolved."
    )
    if outdated:
        print("")
        print(
            f"{INFO}INFO:{RESET} Installing a newer release does not update existing "
            "outputs. Rerunning a directory built with the old release needs "
            "--allow-database-change, or the affected outputs deleted."
        )
    return 0
