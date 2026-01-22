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

rule checkm2:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    params:
        checkm2_module={CHECKM2_MODULE},
        checkm2_db={CHECKM2_DB},
        outdir=f"{OUTPUT_DIR}/dereplicating/checkm2"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Estimating MAG completeness/contamination with CheckM2..."
    shell:
        """
        module load {params.checkm2_module}
        rm -rf {params.outdir}
        mkdir -p {params.outdir}/input {OUTPUT_DIR}/cataloging/final
        for f in {input.genomes}; do
            if [[ "$f" == *.gz ]]; then
                out="{params.outdir}/input/$(basename "$f" .gz)"
                gunzip -c "$f" > "$out"
            else
                ln -sf "$f" "{params.outdir}/input/$(basename "$f")"
            fi
        done
        checkm2 predict --input {params.outdir}/input --output-directory {params.outdir} --threads {threads} --database_path {params.checkm2_db} --force
        python - <<'PY'
            import pandas as pd
            from pathlib import Path

            report_path = Path("{params.outdir}") / "quality_report.tsv"
            df = pd.read_csv(report_path, sep="\\t")
            lower_cols = {{c.lower(): c for c in df.columns}}
            name_col = lower_cols.get("name") or lower_cols.get("genome") or lower_cols.get("bin")
            comp_col = lower_cols.get("completeness")
            cont_col = lower_cols.get("contamination")
            if not name_col or not comp_col or not cont_col:
                raise ValueError(f"CheckM2 output missing required columns: {{df.columns}}")
            out_df = df[[name_col, comp_col, cont_col]].copy()
            out_df.columns = ["genome", "completeness", "contamination"]
            out_df.to_csv("{output.metadata}", index=False)
            PY
        """

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
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
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
        if [ -f "{input.metadata}" ]; then
            # Using existing completeness information
            dRep dereplicate {params.outdir} -p {threads} -g "${{files[@]}}" -sa {params.ani} --genomeInfo {input.metadata}
        else
            # Generate completeness information
            dRep dereplicate {params.outdir} -p {threads} -g "${{files[@]}}" -sa {params.ani}
        fi

        # rename headers in every .fa under dereplicated_genomes/
        for f in {params.outdir}/dereplicated_genomes/*.fa; do
            genome=$(basename "$f" .fa)
            awk -v g="$genome" 'BEGIN{{i=0}} /^>/ {{print ">" g "^" i++; next}} {{print}}' "$f" > tmp && mv tmp "$f"
        done
        """

# Functions to define the input files dynamically.
def get_mag_ids_from_drep(csv_path):
    df = pd.read_csv(csv_path)
    return df["genome"].unique()

def get_derep_final_files(wildcards):
    checkpoint_output = checkpoints.dereplicate.get(**wildcards).output[0]
    selected_bins = get_mag_ids_from_drep(checkpoint_output)
    return expand(f"{OUTPUT_DIR}/dereplicating/final/{{bin_id}}", bin_id=selected_bins)

rule finalize_derep:
    input:
        derep=f"{OUTPUT_DIR}/dereplicating/drep/dereplicated_genomes/{{bin_id}}"
    output:
        f"{OUTPUT_DIR}/dereplicating/final/{{bin_id}}"
    localrule: True
    message: "Collecting dereplicated genomes..."
    shell:
        """
        mkdir -p {OUTPUT_DIR}/dereplicating/final
        cp {input.derep} {output}
        """
