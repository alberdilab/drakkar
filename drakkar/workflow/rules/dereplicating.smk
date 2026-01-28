####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
MASH_MODULE = config["MASH_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
CHECKM2_DB = config["CHECKM2_DB"]
DREP_ANI = float(config.get("DREP_ANI", 0.98))

####
# Workflow rules
####

rule checkm2_report:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        report=f"{OUTPUT_DIR}/dereplicating/checkm2/quality_report.tsv"
    params:
        checkm2_module={CHECKM2_MODULE},
        checkm2_db={CHECKM2_DB},
        outdir=f"{OUTPUT_DIR}/dereplicating/checkm2",
        genome_dir=f"{OUTPUT_DIR}/data/genomes"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Estimating MAG completeness/contamination with CheckM2..."
    shell:
        """
        module load {params.checkm2_module}
        rm -rf {params.outdir}
        rm -rf {params.genome_dir}
        mkdir -p {params.genome_dir}
        for f in {input.genomes}; do
            base=$(basename "$f")
            stem="${{base%.gz}}"
            stem="${{stem%.fa}}"
            stem="${{stem%.fna}}"
            stem="${{stem%.fasta}}"
            out="{params.genome_dir}/${{stem}}.fa"
            if [[ "$f" == *.gz ]]; then
                gunzip -c "$f" > "$out"
            else
                src=$(readlink -f "$f" 2>/dev/null || realpath "$f")
                ln -sf "$src" "$out"
            fi
        done
        checkm2 predict --input {params.genome_dir} --output-directory {params.outdir} --threads {threads} --database_path {params.checkm2_db} --extension fa --force
        rm -rf {params.genome_dir}
        """

rule checkm2_metadata:
    input:
        report=f"{OUTPUT_DIR}/dereplicating/checkm2/quality_report.tsv"
    output:
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    message: "Formatting CheckM2 report for dRep..."
    run:
        import pandas as pd
        import re
        from pathlib import Path

        report_path = Path(input.report)
        df = pd.read_csv(report_path, sep="\t")
        lower_cols = {c.lower(): c for c in df.columns}
        name_col = lower_cols.get("name") or lower_cols.get("genome") or lower_cols.get("bin")
        comp_col = lower_cols.get("completeness")
        cont_col = lower_cols.get("contamination")
        if not name_col or not comp_col or not cont_col:
            raise ValueError(f"CheckM2 output missing required columns: {list(df.columns)}")

        genome_map = {}
        for genome_path in BINS_TO_FILES.values():
            base = Path(genome_path).name
            base_no_gz = re.sub(r"\.gz$", "", base, flags=re.IGNORECASE)
            base_no_ext = re.sub(r"\.(fa|fna|fasta)$", "", base_no_gz, flags=re.IGNORECASE)
            genome_map.setdefault(base_no_ext, base_no_gz)
            genome_map.setdefault(base_no_gz, base_no_gz)

        genomes = []
        missing = []
        for name in df[name_col].astype(str):
            key = name.strip()
            mapped = genome_map.get(key)
            if not mapped:
                missing.append(key)
                genomes.append(key)
            else:
                genomes.append(mapped)

        if missing:
            raise ValueError(f"CheckM2 names not found in data/genomes: {missing}")

        out_df = pd.DataFrame(
            {
                "genome": genomes,
                "completeness": df[comp_col],
                "contamination": df[cont_col],
            }
        )
        Path(output.metadata).parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(output.metadata, index=False)

checkpoint dereplicate:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values()),
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    output:
        Wdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Wdb.csv",
        Cdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Cdb.csv"
    params:
        drep_module={DREP_MODULE},
        mash_module={MASH_MODULE},
        outdir=f"{OUTPUT_DIR}/dereplicating/drep/",
        ani={DREP_ANI},
        uncompressed_dir=f"{OUTPUT_DIR}/dereplicating/drep/uncompressed_genomes"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(60, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.mash_module} {params.drep_module}
        rm -rf {params.outdir}
        mkdir -p {params.uncompressed_dir}
        files=()
        for f in {input.genomes}; do
            if [[ "$f" == *.gz ]]; then
                out="{params.uncompressed_dir}/$(basename "$f" .gz)"
                gunzip -c "$f" > "$out"
                files+=("$out")
            else
                files+=("$f")
            fi
        done
        dRep dereplicate {params.outdir} -p {threads} -g "${{files[@]}}" -sa {params.ani} --genomeInfo {input.metadata}
        """

# Normalize headers in dereplicated genomes with .fa/.fna/.fasta
rule rename_derep_headers:
    input:
        files=lambda wildcards: get_derep_files(wildcards)
    output:
        touch(f"{OUTPUT_DIR}/dereplicating/drep/dereplicated_genomes/.headers_renamed")
    localrule: True
    message: "Renaming headers in dereplicated genomes..."
    shell:
        """
        shopt -s nullglob
        for f in {input.files}; do
            base=$(basename "$f")
            genome="${{base%.fa}}"
            genome="${{genome%.fna}}"
            genome="${{genome%.fasta}}"
            awk -v g="$genome" 'BEGIN{{i=0}} /^>/ {{print ">" g "^" i++; next}} {{print}}' "$f" > tmp && mv tmp "$f"
        done
        touch {output}
        """
# Functions to define the input files dynamically.
def get_mag_ids_from_drep(csv_path):
    df = pd.read_csv(csv_path)
    return df["genome"].unique()

def get_derep_final_files(wildcards):
    checkpoint_output = checkpoints.dereplicate.get(**wildcards).output[0]
    selected_bins = get_mag_ids_from_drep(checkpoint_output)
    return expand(f"{OUTPUT_DIR}/dereplicating/final/{{bin_id}}", bin_id=selected_bins)

def get_derep_files(wildcards):
    checkpoint_output = checkpoints.dereplicate.get(**wildcards).output[0]
    selected_bins = get_mag_ids_from_drep(checkpoint_output)
    return expand(f"{OUTPUT_DIR}/dereplicating/drep/dereplicated_genomes/{{bin_id}}", bin_id=selected_bins)

rule finalize_derep:
    input:
        derep=f"{OUTPUT_DIR}/dereplicating/drep/dereplicated_genomes/{{bin_id}}",
        headers=f"{OUTPUT_DIR}/dereplicating/drep/dereplicated_genomes/.headers_renamed"
    output:
        f"{OUTPUT_DIR}/dereplicating/final/{{bin_id}}"
    localrule: True
    message: "Collecting dereplicated genomes..."
    shell:
        """
        mkdir -p {OUTPUT_DIR}/dereplicating/final
        cp {input.derep} {output}
        """
