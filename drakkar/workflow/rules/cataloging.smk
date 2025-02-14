####
# Define config variables
####

MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

####
# Calculate file sizes
####

preprocess_mb = calculate_file_sizes(PREPROCESS_DIR)
preprocess_mb = {key.replace('_1.fq.gz', ''): value for key, value in preprocess_mb.items()}
preprocess_mb_total = sum(preprocess_mb.values())

####
# Run rules
####

if "individual" in CATALOGING_MODE:
    rule individual_assembly:
        input:
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 16) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 100) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.megahit_module}
            rm -rf {params.outputdir}
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outputdir}
            mv {params.outputdir}/final.contigs.fa {output}
            """

    rule individual_assembly_index:
        input:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna"
        output:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.rev.2.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}"
        threads: 1
        resources:
            mem_mb=32*1024,
            runtime=60
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            """

    rule individual_assembly_map:
        input:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.rev.2.bt2",
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.bam",
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
            """

if "all" in CATALOGING_MODE:
    rule all_assembly:
        input:
            r1=expand(f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz", sample=samples),
            r2=expand(f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz", sample=samples)
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/all"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(1,int(preprocess_mb_total * 1024 / 4 * 2 ** (attempt - 1))),
            runtime=lambda wildcards, attempt: max(15,int((preprocess_mb_total / 4 * 2 ** (attempt - 1))))
        shell:
            """
            module load {params.megahit_module}
            rm -rf {params.outputdir}
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outputdir}
            mv {params.outputdir}/final.contigs.fa {output}
            """

    rule all_assembly_index:
        input:
            f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna"
        output:
            index=f"{OUTPUT_DIR}/cataloging/megahit/all/all.rev.2.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/all"
        threads: 1
        resources:
            mem_mb=32*1024,
            runtime=60
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            """

    rule all_assembly_map:
        input:
            index=f"{OUTPUT_DIR}/cataloging/megahit/all/all.rev.2.bt2",
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/all/{{sample}}.bam",
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/bowtie2/all"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
            """
