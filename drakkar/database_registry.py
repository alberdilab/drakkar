from __future__ import annotations

from pathlib import Path


MANAGED_DATABASES = {
    "kegg": {
        "aliases": ["kofams"],
        "config_key": "KEGG_DB",
        "basename": "kofams",
        "version_label": "current",
        "sources": [
            "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz",
            "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=",
        ],
    },
    "cazy": {
        "aliases": [],
        "config_key": "CAZY_DB",
        "basename": "cazy",
        "version_label": "requested version",
        "sources": [
            "https://bcb.unl.edu/dbCAN2/download/Databases/{version}/dbCAN-HMMdb-{version}.txt",
        ],
    },
    "pfam": {
        "aliases": [],
        "config_key": "PFAM_DB",
        "basename": "pfam",
        "version_label": "current",
        "sources": [
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
            "https://ecdm.loria.fr/data/EC-Pfam_calculated_associations_Extended.csv",
        ],
    },
    "vfdb": {
        "aliases": [],
        "config_key": "VFDB_DB",
        "basename": "vfdb",
        "version_label": "VFDB_setB",
        "sources": [
            "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz",
        ],
    },
    "amr": {
        "aliases": [],
        "config_key": "AMR_DB",
        "basename": "amr",
        "version_label": "NCBIfam-AMRFinder latest",
        "sources": [
            "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz",
            "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv",
        ],
    },
}

DATABASE_ALIASES = {
    alias: name
    for name, definition in MANAGED_DATABASES.items()
    for alias in definition["aliases"]
}


def normalize_managed_database_name(name: str) -> str | None:
    normalized = (name or "").strip().lower()
    if not normalized:
        return None
    normalized = DATABASE_ALIASES.get(normalized, normalized)
    if normalized in MANAGED_DATABASES:
        return normalized
    return None


def database_release_dir(database_name: str, base_directory: str | Path, version: str) -> Path:
    return Path(base_directory) / version


def database_target_path(database_name: str, base_directory: str | Path, version: str) -> Path:
    definition = MANAGED_DATABASES[database_name]
    return database_release_dir(database_name, base_directory, version) / definition["basename"]


def database_sources(database_name: str, version: str | None = None) -> list[str]:
    definition = MANAGED_DATABASES[database_name]
    return [source.format(version=version) for source in definition["sources"]]


def database_source_version_label(database_name: str, version: str | None = None) -> str:
    if database_name == "cazy" and version:
        return f"dbCAN-HMMdb-{version}"
    return MANAGED_DATABASES[database_name]["version_label"]
