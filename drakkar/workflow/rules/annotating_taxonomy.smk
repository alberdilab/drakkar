####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
GTDBTK_MODULE = config["GTDBTK_MODULE"]
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]

# Annotation databases
GTDB_VERSION = str(config.get("gtdb_version", "") or "").strip()
if GTDB_VERSION:
    GTDB_DB_KEY = f"GTDB_DB_{GTDB_VERSION}"
    if GTDB_DB_KEY not in config:
        available = sorted(
            key.replace("GTDB_DB_", "")
            for key in config
            if key.startswith("GTDB_DB_")
        )
        raise ValueError(
            f"GTDB version {GTDB_VERSION} is not configured. "
            f"Add {GTDB_DB_KEY} to config.yaml or choose one of: {', '.join(available)}"
        )
    GTDB_DB = config[GTDB_DB_KEY]
else:
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
        mem_mb = cap_mem_mb(1*1024),
        runtime = cap_runtime(5)
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
        mem_mb=lambda wildcards, attempt: cap_mem_mb(128*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: cap_runtime(120 * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: cap_mem_mb(1*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: cap_runtime(5 * 2 ** (attempt - 1))
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

rule gtdbtk_pruned_trees:
    input:
        batchfile=f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv",
        bacteria_summary=f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    output:
        bacteria=f"{OUTPUT_DIR}/annotating/bacteria.tree",
        ready=touch(f"{OUTPUT_DIR}/annotating/gtdbtk_pruned_trees.done")
    params:
        package_dir={PACKAGE_DIR},
        classify_dir=f"{OUTPUT_DIR}/annotating/gtdbtk/classify",
        archaea_summary=f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        archaea_output=f"{OUTPUT_DIR}/annotating/archaea.tree"
    localrule: True
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: cap_mem_mb(1*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: cap_runtime(5 * 2 ** (attempt - 1))
    message: "Pruning GTDB-Tk trees to keep only input genomes..."
    shell:
        """
        bacteria_tree=""
        for candidate in \
            {params.classify_dir}/gtdbtk.backbone.bac120.classify.tree \
            {params.classify_dir}/gtdbtk.bac120.classify.tree
        do
            if [ -s "$candidate" ]; then
                bacteria_tree="$candidate"
                break
            fi
        done

        if [ -n "$bacteria_tree" ]; then
            python {params.package_dir}/workflow/scripts/prune_gtdbtk_tree.py \
                --input-tree "$bacteria_tree" \
                --batchfile {input.batchfile} \
                --output-tree {output.bacteria}
        else
            : > {output.bacteria}
        fi

        rm -f {params.archaea_output}
        archaea_tree=""
        for candidate in \
            {params.classify_dir}/gtdbtk.backbone.ar53.classify.tree \
            {params.classify_dir}/gtdbtk.ar53.classify.tree
        do
            if [ -s "$candidate" ]; then
                archaea_tree="$candidate"
                break
            fi
        done

        if [ -s {params.archaea_summary} ] && [ -n "$archaea_tree" ]; then
            python {params.package_dir}/workflow/scripts/prune_gtdbtk_tree.py \
                --input-tree "$archaea_tree" \
                --batchfile {input.batchfile} \
                --output-tree {params.archaea_output}
        fi

        touch {output.ready}
        """
