####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
GTDBTK_MODULE = config["GTDBTK_MODULE"]

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
        gtdbtk_module={GTDBTK_MODULE},
        db={GTDB_DB},
        outdir=f"{OUTPUT_DIR}/profiling_genomes/gtdbtk/",
        tmpdir=f"{OUTPUT_DIR}/profiling_genomes/tmp/"
    threads: 24
    resources:
        mem_mb=lambda wildcards, input, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, input, attempt: 120 * 2 ** (attempt - 1)
    message: "Annotating taxonomy using GTDBTK..."
    shell:
        """
        module load {params.gtdbtk_module}
        export GTDBTK_DATA_PATH={params.db}
        mkdir -p {params.tmpdir}
        gtdbtk classify_wf \
            --batchfile {input} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --skip_ani_screen
        """
