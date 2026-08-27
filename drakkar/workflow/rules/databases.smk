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
FOLDSEEK_MODULE = config["FOLDSEEK_MODULE"]

DATABASE_NAME = normalize_managed_database_name(config.get("database_name", ""))
DATABASE_DIRECTORY = config.get("database_directory", "")
DATABASE_VERSION = config.get("database_version", "")
DOWNLOAD_RUNTIME = int(config.get("database_download_runtime", 120))
INSTALL_DIR = Path(OUTPUT_DIR)

if not DATABASE_NAME or DATABASE_NAME not in MANAGED_DATABASES:
    raise ValueError(f"Unsupported database_name: {config.get('database_name')}")
if DOWNLOAD_RUNTIME <= 0:
    raise ValueError(f"database_download_runtime must be positive, got {DOWNLOAD_RUNTIME}")

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
            ko_list=f"{TARGET_DB}_ko_list.tsv",
            archive_url=DATABASE_SOURCES[0],
            json_url=DATABASE_SOURCES[1],
            ko_list_url=DATABASE_SOURCES[2],
            hmmer_module=HMMER_MODULE
        threads: 1
        resources:
            runtime=lambda wildcards, attempt: cap_runtime(DOWNLOAD_RUNTIME * 2 ** (attempt - 1))
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.json}" "{params.ko_list}" "{params.ko_list}.gz" "{params.archive}"
            rm -rf "{OUTPUT_DIR}/profiles"
            if ! curl -L --fail --output "{params.archive}" "{params.archive_url}"; then
                echo "Unable to download KEGG KOfam archive: {params.archive_url}" >&2
                echo "If this archive version does not exist, check available versions at https://www.genome.jp/ftp/db/kofam/archives/" >&2
                rm -f "{params.archive}"
                exit 1
            fi
            curl -L --fail --output "{params.json}" "{params.json_url}"
            curl -L --fail --output "{params.ko_list}.gz" "{params.ko_list_url}"
            gunzip -f "{params.ko_list}.gz"
            tar -xzf "{params.archive}" -C "{OUTPUT_DIR}"
            find "{OUTPUT_DIR}/profiles" -type f -name "*.hmm" | sort | xargs cat > "{params.db}"
            rm -f "{params.archive}"
            rm -rf "{OUTPUT_DIR}/profiles"
            module purge
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
            if ! curl -L --fail --output "{params.db}.tmp" "{params.url}"; then
                echo "Unable to download CAZy dbCAN database: {params.url}" >&2
                echo "If this dbCAN release does not exist, check available versions at https://pro.unl.edu/dbCAN2/browse_download.php" >&2
                rm -f "{params.db}.tmp"
                exit 1
            fi
            mv "{params.db}.tmp" "{params.db}"
            module purge
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
            hmm_url=DATABASE_SOURCES[0],
            ec_url=DATABASE_SOURCES[1],
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.gz" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.ec}"
            if ! curl -L --fail --output "{params.db}.gz" "{params.hmm_url}"; then
                echo "Unable to download Pfam release: {params.hmm_url}" >&2
                echo "If this Pfam release does not exist, check available versions at https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" >&2
                rm -f "{params.db}.gz"
                exit 1
            fi
            gunzip -f "{params.db}.gz"
            curl -L --fail --output "{params.ec}" "{params.ec_url}"
            module purge
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
            module purge
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
            archive_url=DATABASE_SOURCES[0],
            tsv_url=DATABASE_SOURCES[1],
            hmmer_module=HMMER_MODULE
        threads: 1
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -f "{params.db}" "{params.db}.h3f" "{params.db}.h3i" "{params.db}.h3m" "{params.db}.h3p" "{params.archive}" "{params.tsv}"
            rm -rf "{params.hmmdir}"
            if ! curl -L --fail --output "{params.archive}" "{params.archive_url}"; then
                echo "Unable to download AMRFinder release: {params.archive_url}" >&2
                echo "If this AMRFinder release does not exist, check available versions at https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/" >&2
                rm -f "{params.archive}"
                exit 1
            fi
            tar -xzf "{params.archive}" -C "{OUTPUT_DIR}"
            find "{params.hmmdir}" -type f -name "*.HMM" | sort | xargs cat > "{params.db}"
            if ! curl -L --fail --output "{params.tsv}" "{params.tsv_url}"; then
                echo "Unable to download AMRFinder metadata table: {params.tsv_url}" >&2
                echo "If this AMRFinder release does not exist, check available versions at https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/" >&2
                rm -f "{params.tsv}" "{params.archive}"
                exit 1
            fi
            module purge
            module load {params.hmmer_module}
            hmmpress -f "{params.db}"
            rm -f "{params.archive}"
            rm -rf "{params.hmmdir}"
            touch {output}
            """

if DATABASE_NAME == "amrfinderplus":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/amrfinderplus.done")
        params:
            source=DATABASE_SOURCES[0],
            staging=f"{OUTPUT_DIR}/amrfinderplus.download",
            expected_version=DATABASE_VERSION
        threads: 1
        conda:
            f"{PACKAGE_DIR}/workflow/envs/amr_amrfinder.yaml"
        resources:
            runtime=lambda wildcards, attempt: cap_runtime(DOWNLOAD_RUNTIME * 2 ** (attempt - 1))
        shell:
            r"""
            set -euo pipefail
            rm -rf {params.staging:q}
            mkdir -p {params.staging:q}

            for filename in \
                AMR.LIB \
                AMRProt.fa \
                AMRProt-mutation.tsv \
                AMRProt-suppress.tsv \
                AMRProt-susceptible.fa \
                AMRProt-susceptible.tsv \
                AMR_CDS.fa \
                database_format_version.txt \
                fam.tsv \
                taxgroup.tsv \
                version.txt \
                changes.txt
            do
                curl -L --fail \
                    --output "{params.staging}/$filename" \
                    "{params.source}$filename"
            done

            while read -r taxgroup gpipe mutation_count rest
            do
                case "$taxgroup" in
                    ""|\#*) continue ;;
                esac
                if [ "$mutation_count" -gt 0 ] 2>/dev/null
                then
                    for extension in fa tsv
                    do
                        filename="AMR_DNA-${{taxgroup}}.${{extension}}"
                        curl -L --fail \
                            --output "{params.staging}/$filename" \
                            "{params.source}$filename"
                    done
                fi
            done < {params.staging:q}/taxgroup.tsv

            if ! grep -Fxq {params.expected_version:q} {params.staging:q}/version.txt
            then
                echo "Downloaded AMRFinderPlus database version does not match requested version {params.expected_version}." >&2
                cat {params.staging:q}/version.txt >&2
                exit 1
            fi

            amrfinder_index {params.staging:q}
            test -s {params.staging:q}/AMRProt.fa.phr
            test -s {params.staging:q}/AMRProt.fa.pin
            test -s {params.staging:q}/AMRProt.fa.psq

            cp -a {params.staging:q}/. {OUTPUT_DIR:q}/
            rm -rf {params.staging:q}
            touch {output}
            """

if DATABASE_NAME == "card":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/card.done")
        params:
            archive=f"{OUTPUT_DIR}/card-data.tar.bz2",
            staging=f"{OUTPUT_DIR}/card.download",
            url=DATABASE_SOURCES[0],
            validator=f"{PACKAGE_DIR}/workflow/scripts/validate_card_database.py",
            expected_version=DATABASE_VERSION
        threads: 1
        conda:
            f"{PACKAGE_DIR}/workflow/envs/amr_rgi.yaml"
        resources:
            runtime=lambda wildcards, attempt: cap_runtime(DOWNLOAD_RUNTIME * 2 ** (attempt - 1))
        shell:
            r"""
            set -euo pipefail
            rm -f {params.archive:q}
            rm -rf {params.staging:q} {OUTPUT_DIR:q}/localDB
            mkdir -p {params.staging:q}

            curl -L --fail --output {params.archive:q} {params.url:q}
            tar -xjf {params.archive:q} -C {params.staging:q} ./card.json
            python {params.validator:q} \
                {params.staging:q}/card.json \
                --expect {params.expected_version:q}
            mv {params.staging:q}/card.json {TARGET_DB:q}
            rm -rf {params.staging:q}
            rm -f {params.archive:q}

            cd {OUTPUT_DIR:q}
            rgi load --card_json {TARGET_DB:q} --local
            test -s {OUTPUT_DIR:q}/localDB/card.json
            test -s {OUTPUT_DIR:q}/localDB/loaded_databases.json
            touch {output}
            """

if DATABASE_NAME == "foldseek":
    rule prepare_database:
        output:
            touch(f"{OUTPUT_DIR}/foldseek.done")
        params:
            foldseek_db=str(TARGET_DB),
            prostt5=f"{OUTPUT_DIR}/prostt5",
            sprot=f"{OUTPUT_DIR}/uniprot_sprot.dat.gz",
            mapping=f"{OUTPUT_DIR}/foldseek_map.tsv",
            tmp=f"{OUTPUT_DIR}/tmp",
            sprot_url=DATABASE_SOURCES[2],
            foldseek_module=FOLDSEEK_MODULE,
            package_dir=PACKAGE_DIR
        threads: 1
        conda:
            f"{PACKAGE_DIR}/workflow/envs/annotating_function.yaml"
        resources:
            runtime=lambda wildcards, attempt: cap_runtime(DOWNLOAD_RUNTIME * 2 ** (attempt - 1))
        shell:
            """
            set -euo pipefail
            mkdir -p "{OUTPUT_DIR}"
            rm -rf "{params.tmp}" "{params.prostt5}"*
            find "{OUTPUT_DIR}" -maxdepth 1 -type f -name "foldseek_db*" -delete
            rm -f "{params.mapping}" "{params.sprot}"
            module load {params.foldseek_module}
            # Pre-built AlphaFold/Swiss-Prot structure DB and the ProstT5 weights
            # are fetched and formatted by foldseek itself from its mirror.
            foldseek databases Alphafold/Swiss-Prot "{params.foldseek_db}" "{params.tmp}"
            foldseek databases ProstT5 "{params.prostt5}" "{params.tmp}"
            # UniProt accession -> KO/EC/Pfam map for interpreting structural hits.
            curl -L --fail --output "{params.sprot}" "{params.sprot_url}"
            PYTHON_BIN="${{CONDA_PREFIX}}/bin/python"
            "$PYTHON_BIN" "{params.package_dir}/workflow/scripts/build_foldseek_function_map.py" \
                -i "{params.sprot}" \
                -o "{params.mapping}"
            rm -f "{params.sprot}"
            rm -rf "{params.tmp}"
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
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}.json"), Path(f"{TARGET_DB}_ko_list.tsv")]
        elif DATABASE_NAME == "cazy":
            checksums = [Path(str(TARGET_DB))]
        elif DATABASE_NAME == "pfam":
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}_ec.tsv")]
        elif DATABASE_NAME == "vfdb":
            checksums = [Path(f"{TARGET_DB}.idx"), Path(f"{TARGET_DB}.tsv")]
        elif DATABASE_NAME == "amr":
            checksums = [Path(str(TARGET_DB)), Path(f"{TARGET_DB}.tsv")]
        elif DATABASE_NAME == "amrfinderplus":
            checksums = [
                Path(str(TARGET_DB)),
                Path(f"{OUTPUT_DIR}/version.txt"),
                Path(f"{OUTPUT_DIR}/database_format_version.txt"),
                Path(f"{OUTPUT_DIR}/fam.tsv"),
            ]
        elif DATABASE_NAME == "card":
            checksums = [
                Path(str(TARGET_DB)),
                Path(f"{OUTPUT_DIR}/localDB/card.json"),
                Path(f"{OUTPUT_DIR}/localDB/loaded_databases.json"),
            ]
        elif DATABASE_NAME == "foldseek":
            checksums = [Path(str(TARGET_DB)), Path(f"{OUTPUT_DIR}/foldseek_map.tsv")]

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
