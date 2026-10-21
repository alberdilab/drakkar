####
# Per-sample subsampling for standalone cataloging and expressing.
# Only multi-lane samples need a temporary concatenation before seqkit.
####

rule combine_downstream_reads:
    wildcard_constraints:
        mate="[12]"
    input:
        lambda wildcards: (DOWNSTREAM_READS1 if wildcards.mate == "1" else DOWNSTREAM_READS2)[wildcards.sample]
    output:
        temp(f"{OUTPUT_DIR}/{WORKFLOW}/subsample_input/{{sample}}_{{mate}}.fq.gz")
    log:
        f"{OUTPUT_DIR}/{WORKFLOW}/subsample_input/{{sample}}_{{mate}}.log"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 256)) * 2 ** (attempt - 1))
    shell:
        """
        cat {input:q} > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """

rule subsample_downstream_reads:
    input:
        r1=lambda wildcards: downstream_subsample_input(wildcards, 1),
        r2=lambda wildcards: downstream_subsample_input(wildcards, 2)
    output:
        r1=f"{OUTPUT_DIR}/{WORKFLOW}/subsample/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/{WORKFLOW}/subsample/{{sample}}_2.fq.gz",
        stats=f"{OUTPUT_DIR}/{WORKFLOW}/subsample/{{sample}}_subsample.tsv"
    log:
        f"{OUTPUT_DIR}/{WORKFLOW}/subsample/{{sample}}.log"
    params:
        seqkit_module=config["SEQKIT_MODULE"],
        script=f"{PACKAGE_DIR}/workflow/scripts/subsample_reads.py",
        target=lambda wildcards: SAMPLE_TO_SUBSAMPLE[wildcards.sample],
        seed=11
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1024, int(input.size_mb / 8)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 256)) * 2 ** (attempt - 1))
    message: "Subsampling {WORKFLOW} read pairs for {wildcards.sample}..."
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
