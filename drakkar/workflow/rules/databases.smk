####
# Config variables
####

from datetime import datetime, timezone
import importlib.util
from pathlib import Path

PACKAGE_DIR = config["package_dir"]
registry_path = Path(PACKAGE_DIR) / "database_registry.py"
registry_spec = importlib.util.spec_from_file_location("drakkar_database_registry", registry_path)
if registry_spec is None or registry_spec.loader is None:
    raise ImportError(f"Unable to load database registry from {registry_path}")

database_registry = importlib.util.module_from_spec(registry_spec)
registry_spec.loader.exec_module(database_registry)

MANAGED_DATABASES = database_registry.MANAGED_DATABASES
database_target_path = database_registry.database_target_path
database_sources = database_registry.database_sources
database_source_version_label = database_registry.database_source_version_label
normalize_managed_database_name = database_registry.normalize_managed_database_name

HMMER_MODULE = config["HMMER_MODULE"]
MMSEQS2_MODULE = config["MMSEQS2_MODULE"]

DATABASE_NAME = normalize_managed_database_name(config.get("database_name", ""))
DATABASE_DIRECTORY = config.get("database_directory", "")
DATABASE_VERSION = config.get("database_version", "")
INSTALL_DIR = Path(OUTPUT_DIR)

if not DATABASE_NAME or DATABASE_NAME not in MANAGED_DATABASES:
    raise ValueError(f"Unsupported database_name: {config.get('database_name')}")

DATABASE_DEFINITION = MANAGED_DATABASES[DATABASE_NAME]
TARGET_DB = database_target_path(DATABASE_NAME, DATABASE_DIRECTORY, DATABASE_VERSION)
DATABASE_SOURCES = database_sources(DATABASE_NAME, DATABASE_VERSION)
SOURCE_VERSION = database_source_version_label(DATABASE_NAME, DATABASE_VERSION)

if INSTALL_DIR.resolve() != TARGET_DB.parent.resolve():
    raise ValueError(
        f"Database install directory mismatch: output_dir={INSTALL_DIR} but target release dir is {TARGET_DB.parent}"
    )

####
# Database preparation rules
####

if DATABASE_NAME == "kegg":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/kegg.done")
        params:
            db=str(TARGET_DB),
            archive=f"{OUTPUT_DIR}/profiles.tar.gz",
            json=f"{TARGET_DB}.json",
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.json}" "{params.archive}"
            rm -rf "{OUTPUT_DIR}/profiles"
            wget -O "{params.archive}" "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
            curl -L -o "{params.json}" "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir="
            tar -xzf "{params.archive}" -C "{OUTPUT_DIR}"
            find "{OUTPUT_DIR}/profiles" -type f -name "*.hmm" | sort | xargs cat > "{params.db}"
            rm -rf "{OUTPUT_DIR}/profiles"
            module load {params.hmmer_module}
            hmmpress -f "{params.db}"
            touch {output}
            """

if DATABASE_NAME == "cazy":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/cazy.done")
        params:
            db=str(TARGET_DB),
            url=DATABASE_SOURCES[0],
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.tmp" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p"
            curl -L --fail --output "{params.db}.tmp" "{params.url}"
            mv "{params.db}.tmp" "{params.db}"
            module load {params.hmmer_module}
            hmmpress -f "{params.db}"
            touch {output}
            """

if DATABASE_NAME == "pfam":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/pfam.done")
        params:
            db=str(TARGET_DB),
            ec=f"{TARGET_DB}_ec.tsv",
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.gz" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.ec}"
            wget -O "{params.db}.gz" "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
            gunzip -f "{params.db}.gz"
            wget -O "{params.ec}" "https://ecdm.loria.fr/data/EC-Pfam_calculated_associations_Extended.csv"
            module load {params.hmmer_module}
            hmmpress -f "{params.db}"
            touch {output}
            """

if DATABASE_NAME == "vfdb":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/vfdb.done")
        params:
            db=str(TARGET_DB),
            db_prefix=TARGET_DB.name,
            fasta=f"{OUTPUT_DIR}/VFDB_setB_pro.fas",
            mapping=f"{TARGET_DB}.tsv",
            tmp=f"{OUTPUT_DIR}/tmp",
            mmseqs2_module=MMSEQS2_MODULE,
            package_dir=PACKAGE_DIR
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.fasta}" "{params.fasta}.gz" "{params.mapping}"
            rm -rf "{params.tmp}"
            find "{OUTPUT_DIR}" -maxdepth 1 -type f -name "{params.db_prefix}*" -delete
            wget -O "{params.fasta}.gz" "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"
            gunzip -f "{params.fasta}.gz"
            python "{params.package_dir}/workflow/scripts/get_fvid.py" "{params.fasta}" "{params.mapping}"
            module load {params.mmseqs2_module}
            mmseqs createdb "{params.fasta}" "{params.db}"
            mmseqs createindex "{params.db}" "{params.tmp}"
            rm -rf "{params.tmp}"
            touch {output}
            """

if DATABASE_NAME == "amr":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/amr.done")
        params:
            db=str(TARGET_DB),
            archive=f"{OUTPUT_DIR}/NCBIfam-AMRFinder.HMM.tar.gz",
            tsv=f"{TARGET_DB}.tsv",
            hmmdir=f"{OUTPUT_DIR}/HMM",
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.archive}" "{params.tsv}"
            rm -rf "{params.hmmdir}"
            wget -O "{params.archive}" "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz"
            tar -xzf "{params.archive}" -C "{OUTPUT_DIR}"
            find "{params.hmmdir}" -type f -name "*.HMM" | sort | xargs cat > "{params.db}"
            wget -O "{params.tsv}" "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv"
            module load {params.hmmer_module}
            hmmpress -f "{params.db}"
            rm -rf "{params.hmmdir}"
            touch {output}
            """

rule write_database_versions:
    input:
        f"{OUTPUT_DIR}/{DATABASE_NAME}.done"
    output:
        f"{OUTPUT_DIR}/database_versions.yaml"
    run:
        import hashlib
        import yaml

        def sha256sum(path):
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            return digest.hexdigest()

        checksums = []
        if DATABASE_NAME == "kegg":
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}.json")]
        elif DATABASE_NAME == "cazy":
            checksums = [Path(str(TARGET_DB))]
        elif DATABASE_NAME == "pfam":
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}_ec.tsv")]
        elif DATABASE_NAME == "vfdb":
            checksums = [Path(f"{TARGET_DB}.idx"), Path(f"{TARGET_DB}.tsv")]
        elif DATABASE_NAME == "amr":
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}.tsv")]

        file_info = []
        for target in checksums:
            if target.exists():
                file_info.append(
                    {
                        "path": str(target),
                        "sha256": sha256sum(target),
                        "size_bytes": target.stat().st_size,
                    }
                )

        version_info = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "database": DATABASE_NAME,
            "config_key": DATABASE_DEFINITION["config_key"],
            "base_directory": DATABASE_DIRECTORY,
            "release_directory": str(INSTALL_DIR),
            "requested_version": DATABASE_VERSION,
            "source_version": SOURCE_VERSION,
            "sources": DATABASE_SOURCES,
            "default_target": str(TARGET_DB),
            "files": file_info,
        }

        with open(output[0], "w") as handle:
            yaml.safe_dump(version_info, handle, sort_keys=False)
