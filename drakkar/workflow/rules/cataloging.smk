####
# Define config variables
####

MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

####
# Calculate file sizes
####

reads_mb = calculate_file_sizes(READS_DIR)
reads_mb = {key.replace('_1.fq.gz', ''): value for key, value in reads_mb.items()}
reads_mb_total = sum(reads_mb.values())

####
# Run rules
####

if CATALOGING_MODE == "individual":
    rule individual_assembly:
        input:
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}"
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

    rule assembly_index:
        input:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}.fna"
        output:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}.rev.2.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}"
        threads: 1
        resources:
            mem_mb=32*1024,
            runtime=60
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            """

    rule assembly_map:
        input:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}.rev.2.bt2",
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}.bam",
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
            """

if CATALOGING_MODE == "coassembly":
    rule coassembly:
        input:
            r1=expand(f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz", sample=samples),
            r2=expand(f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz", sample=samples)
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/coassembly.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/"
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
