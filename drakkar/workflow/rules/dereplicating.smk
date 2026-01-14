####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
MASH_MODULE = config["MASH_MODULE"]
DREP_ANI = float(config.get("DREP_ANI", 0.98))

####
# Workflow rules
####

checkpoint dereplicate:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        Wdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Wdb.csv",
        Cdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Cdb.csv"
    params:
        drep_module={DREP_MODULE},
        mash_module={MASH_MODULE},
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv",
        outdir=f"{OUTPUT_DIR}/dereplicating/drep/",
        ani={DREP_ANI}
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.mash_module} {params.drep_module}
        rm -rf {params.outdir}
        if [ -f "{params.metadata}" ]; then
            # Using existing completeness information
            dRep dereplicate {params.outdir} -p {threads} -g {input.genomes} -sa {params.ani} --genomeInfo {params.metadata}
        else
            # Generate completeness information
            dRep dereplicate {params.outdir} -p {threads} -g {input.genomes} -sa {params.ani}
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
