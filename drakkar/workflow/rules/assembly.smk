####
# Define config variables
####

MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
ASSEMBLY_MODE = config["ASSEMBLY_MODE"]

if ASSEMBLY_MODE == "individual":
    rule individual_assembly:
        input:
            r1=f"{READS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{READS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/assembly/final/{{sample}}.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/assembly/megahit/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 8) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.megahit_module}
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outputdir}
            mv {params.outputdir}/final.contigs.fa {output}
            """

if ASSEMBLY_MODE == "coassembly":
    rule coassembly:
        input:
            r1=expand(f"{READS_DIR}/{{sample}}_1.fq.gz", sample=samples)
            r2=expand(f"{READS_DIR}/{{sample}}_2.fq.gz", sample=samples)
        output:
            f"{OUTPUT_DIR}/assembly/final/coassembly.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/assembly/megahit/"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 8) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.megahit_module}
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outputdir}
            mv {params.outputdir}/final.contigs.fa {output}
            """
