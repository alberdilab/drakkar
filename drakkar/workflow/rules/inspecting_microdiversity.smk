####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
COVERM_MODULE = config["COVERM_MODULE"]
GATK_MODULE = config["GATK_MODULE"]
SEQTK_MODULE = config["SEQTK_MODULE"]

####
# Workflow rules
####

rule subset_catalogue:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    output:
        f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fna"
    threads: 1
    params:
        seqtk_module={SEQTK_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.seqtk_module}
        printf "%s\n" "{wildcards.genome}" | seqtk subseq {input} - > {output}
        """

rule subset_bam:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam"
    output:
        bam=f"{OUTPUT_DIR}/inspecting_microdiversity/samples/{{genome}}/{{sample}}.bam",
        bai=f"{OUTPUT_DIR}/inspecting_microdiversity/samples/{{genome}}/{{sample}}.bam.bai"
    threads: 1
    params:
        samtools_module={SAMTOOLS_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.samtools_module}
        samtools view -h {input} {wildcards.genome} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule remove_duplicates:
    input:
        f"{OUTPUT_DIR}/inspecting_microdiversity/samples/{{genome}}/{{sample}}.bam"
    output:
        dedup=temp(f"{OUTPUT_DIR}/inspecting_microdiversity/gatk/{{genome}}/{{sample}}.tmp.bam"),
        clean=f"{OUTPUT_DIR}/inspecting_microdiversity/gatk/{{genome}}/{{sample}}.bam",
        metrics=f"{OUTPUT_DIR}/inspecting_microdiversity/gatk/{{genome}}/{{sample}}.txt"
    params:
        picard_module={PICARD_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        tempdir=f"{OUTPUT_DIR}/inspecting_microdiversity/gatk/{{genome}}/{{sample}}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load python/3.12.8 openjdk/17.0.8 {params.gatk_module}
        gatk MarkDuplicates INPUT={input} OUTPUT={output.dedup} REMOVE_DUPLICATES=True METRICS_FILE={output.metrics}
        gatk CleanSam INPUT={output.dedup} OUTPUT={output.clean} TMP_DIR={params.tempdir}
        """

rule create_dict:
    input:
        f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fna"
    output:
        fai=f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fai",
        dict=f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.dict"
    params:
        picard_module={PICARD_MODULE},
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load python/3.12.8 openjdk/17.0.8 {params.gatk_module} {params.samtools_module}
        samtools faidx {input} -o {output.fai}
        gatk CreateSequenceDictionary R={input} O={output.dict}
        """

rule call_variants:
    input:
        catalogue=f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fna",
        fai=f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fai",
        bam=f"{OUTPUT_DIR}/inspecting_microdiversity/gatk/{{genome}}/{{sample}}.bam"
    output:
        f"{OUTPUT_DIR}/profiling_genomes/gatk/{{genome}}/{{sample}}.vcf.gz"
    params:
        gatk_module={GATK_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Calling SNVs in sample {wildcards.sample}..."
    shell:
        """
        module load python/3.12.8 openjdk/17.0.8 {params.gatk_module}
		gatk HaplotypeCaller -ploidy=1 -R {input.catalogue} -I {input.bam} -O {output}
        """

rule merge_variants:
    input:
        catalogue=f"{OUTPUT_DIR}/inspecting_microdiversity/catalogue/{{genome}}.fna",
        gvcfs=expand(f"{OUTPUT_DIR}/profiling_genomes/gatk/{{genome}}/{{sample}}.vcf.gz", sample=samples)
    output:
        f"{OUTPUT_DIR}/inspecting_microdiversity/variants/{{genome}}.vcf.gz"
    params:
        gatk_module={GATK_MODULE}
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Merging haplotypes across samples..."
    shell:
        """
        module load python/3.12.8 openjdk/17.0.8 {params.gatk_module}
		gatk combinegvcfs {input.gvcfs} -R {input.catalogue} -O {output}
        """
