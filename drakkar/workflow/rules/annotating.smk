####
# Define config variables
####

PACKAGE_DIR = config["package_dir"]
GTDBTK_MODULE = config["GTDBTK_MODULE"]

rule gtdbtk_input:
    input:
        json=f"{OUTPUT_DIR}/data/bins_to_files.json",
        names=expand("{bin_name}", bin_path=BINS_TO_FILES.keys()),
        paths=expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    params:
        package_dir={PACKAGE_DIR}
    localrule: True
    threads: 1
    resources:
        mem_mb = 1*1024,
        runtime = 5
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/gtdbtk_input.py --names {input.names} --paths {input.paths} --output {output}
        """

rule gtdbtk:
    input:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        gtdbtk_module={GTDBTK_MODULE},
        outdir=f"{OUTPUT_DIR}/profiling_genomes/gtdbtk/",
        tmpdir=f"{OUTPUT_DIR}/profiling_genomes/gtdbtk/tmp"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.gtdbtk_module}
        gtdbtk classify_wf \
            --batchfile {input} \
            --out_dir {params.outdir} \
            --tmpdir {params.tmpdir} \
            --cpus {threads} \
            --skip_ani_screen
        """
