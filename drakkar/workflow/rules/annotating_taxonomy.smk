####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
GTDBTK_MODULE = config["GTDBTK_MODULE"]
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]

# Annotation databases
GTDB_DB = config["GTDB_DB"]

####
# Workflow rules
####

rule gtdbtk_input:
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    params:
        package_dir={PACKAGE_DIR},
        names=expand("{mag_name}", mag_name=MAGS_TO_FILES.keys()),
        paths=expand("{mag_path}", mag_path=MAGS_TO_FILES.values())
    localrule: True
    threads: 1
    resources:
        mem_mb = 1*1024,
        runtime = 5
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/gtdbtk_input.py --names {params.names} --paths {params.paths} --output {output}
        """

rule gtdbtk:
    input:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        prodigal_module={PRODIGAL_MODULE},
        gtdbtk_module={GTDBTK_MODULE},
        db={GTDB_DB},
        outdir=f"{OUTPUT_DIR}/annotating/gtdbtk/",
        tmpdir=f"{OUTPUT_DIR}/annotating/tmp/"
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_taxonomy.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 120 * 2 ** (attempt - 1)
    message: "Annotating taxonomy using GTDBTK..."
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        mkdir -p {params.tmpdir}
        gtdbtk classify_wf \
            --batchfile {input} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --skip_ani_screen
        """

rule gtdbtk_table:
    input:
        f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    output:
        f"{OUTPUT_DIR}/annotating/genome_taxonomy.tsv"
    params:
        archaea=f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        archaea2=f"{OUTPUT_DIR}/annotating/gtdbtk/classify/ar53.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 5 * 2 ** (attempt - 1)
    message: "Annotating taxonomy using GTDBTK..."
    shell:
        """
        # Create a merged summary output for DRAM:
        if [ -s {params.archaea} ]
        then
        sed '1d;' {params.archaea} > {params.archaea2}
        cat {input} {params.archaea2} > {output}
        rm {params.archaea2}
        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {input} > {output}
        fi
        """
