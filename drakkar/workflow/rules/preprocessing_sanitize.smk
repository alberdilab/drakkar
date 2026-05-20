####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
SEQKIT_MODULE = config["SEQKIT_MODULE"]

####
# Workflow rules
####

# Sanitize each read file independently with seqkit sana to remove malformed records,
# then re-pair the two files so only reads present in both are retained.
# Outputs sanitized paired FASTQ files and a stats file recording the retained read count.

rule seqkit_sana:
    input:
        r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
    output:
        r1=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_2.fq.gz",
        stats=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_sana.tsv"
    params:
        seqkit_module=SEQKIT_MODULE,
        workdir=lambda wildcards: f"{OUTPUT_DIR}/preprocessing/seqkit/.tmp_{wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 512) * 2 ** (attempt - 1)))
    message: "Sanitizing FASTQ files for sample {wildcards.sample}..."
    shell:
        """
        module purge
        module load {params.seqkit_module}
        mkdir -p {params.workdir}/paired
        seqkit sana {input.r1} -o {params.workdir}/sana_1.fq.gz
        seqkit sana {input.r2} -o {params.workdir}/sana_2.fq.gz
        seqkit pair \
            -1 {params.workdir}/sana_1.fq.gz \
            -2 {params.workdir}/sana_2.fq.gz \
            -O {params.workdir}/paired
        mv {params.workdir}/paired/sana_1.fq.gz {output.r1}
        mv {params.workdir}/paired/sana_2.fq.gz {output.r2}
        count=$(zcat {output.r1} | awk 'NR%4==1{{count++}} END{{print count * 2}}')
        printf "sample\tseqkit_sana_reads\n{wildcards.sample}\t$count" > {output.stats}
        rm -rf {params.workdir}
        """
