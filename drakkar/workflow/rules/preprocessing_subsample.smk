####
# Subsample paired raw reads before sanitizing or quality filtering.
####

rule subsample_reads:
    input:
        r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
    output:
        r1=f"{OUTPUT_DIR}/preprocessing/subsample/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/subsample/{{sample}}_2.fq.gz",
        stats=f"{OUTPUT_DIR}/preprocessing/subsample/{{sample}}_subsample.tsv"
    log:
        f"{OUTPUT_DIR}/preprocessing/subsample/{{sample}}.log"
    params:
        seqkit_module=config["SEQKIT_MODULE"],
        script=f"{PACKAGE_DIR}/workflow/scripts/subsample_reads.py",
        target=SUBSAMPLE,
        seed=11
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1024, int(input.size_mb / 8)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 256)) * 2 ** (attempt - 1))
    message: "Subsampling raw read pairs for {wildcards.sample}..."
    shell:
        """
        py=$(command -v python3 || command -v python)
        module purge
        module load {params.seqkit_module}
        "$py" {params.script:q} \
            --sample {wildcards.sample:q} --target {params.target:q} --seed {params.seed} \
            --r1 {input.r1:q} --r2 {input.r2:q} \
            --out-r1 {output.r1:q} --out-r2 {output.r2:q} \
            --stats {output.stats:q} > {log:q} 2>&1
        """
