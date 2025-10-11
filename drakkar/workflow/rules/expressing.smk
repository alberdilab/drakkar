####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
SUBREAD_MODULE=config["SUBREAD_MODULE"]

####
# Workflow rules
####

rule prodigal:
    input:
        lambda wildcards: MAGS_TO_FILES[wildcards.mag]
    output:
        nt=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.fna",
        aa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa",
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff"
    params:
        prodigal_module={PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    threads: 1
    message: "Predicting genes of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.prodigal_module}
        if [[ "{input}" == *.gz ]]; then
            gzip -dc {input} | prodigal -d {output.nt} -a {output.aa} -o {output.gff} -f gff
        else
            prodigal -i {input} -d {output.nt} -a {output.aa} -o {output.gff} -f gff
        fi
        """

rule merge_metagenome:
    input:
        MAGS_TO_FILES.values()
    output:
        f"{OUTPUT_DIR}/expressing/reference/metagenome.fna"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        cat {input} > {output}
        """

rule merge_metagenome_gff:
    input:
        expand(f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff", mag=mags)
    output:
        f"{OUTPUT_DIR}/expressing/reference/metagenome.gff"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        cat {input} > {output}
        """

rule index_metagenome:
    input:
        f"{OUTPUT_DIR}/expressing/reference/metagenome.fna"
    output:
        index=f"{OUTPUT_DIR}/expressing/reference/metagenome.rev.1.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/expressing/reference/metagenome"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        """

rule map_to_metagenome:
    input:
        index=f"{OUTPUT_DIR}/expressing/reference/metagenome.rev.1.bt2",
        r1=lambda wildcards: TRANSCRIPTOME_TO_READS1[wildcards.sample],
        r2=lambda wildcards: TRANSCRIPTOME_TO_READS2[wildcards.sample]
    output:
        f"{OUTPUT_DIR}/expressing/bowtie2/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=f"{OUTPUT_DIR}/expressing/reference/metagenome"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against genome catalogue..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """
    
rule quantify:
    input:
        bams=expand(f"{OUTPUT_DIR}/expressing/bowtie2/{{sample}}.bam", sample=samples),
        annotation=f"{OUTPUT_DIR}/expressing/reference/metagenome.gff"
    output:
        counts=f"{OUTPUT_DIR}/expressing/featurecounts/counts.tsv",
        summary=f"{OUTPUT_DIR}/expressing/featurecounts/counts.tsv.summary"
    threads: 8
    params:
        subread_module={SUBREAD_MODULE},
        extra="-F GFF -t CDS,tRNA,rRNA -g ID -p",
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 2) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.subread_module}
        featureCounts {params.extra} -a {input.annotation} -o {output.counts} {input.bams}
        """