####
# Variables parsed from the config.yaml file
####

PACKAGE_DIR = config["package_dir"]

# Software modules
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
BCFTOOLS_MODULE = config["BCFTOOLS_MODULE"]
GATK_MODULE = config["GATK_MODULE"]
JAVA_MODULE = config.get("JAVA_MODULE", "openjdk/17.0.8")
PYTHON_MODULE = config.get("PYTHON_MODULE", "python/3.12.8")

# Number of interval shards HaplotypeCaller is scattered over. Host references have few
# but very large chromosomes, so scattering per contig would leave one job doing most of
# the work; SplitIntervals subdivides contigs to produce evenly sized shards instead.
SCATTER_COUNT = int(config.get("GENOMICS_SCATTER_COUNT", 24))
SHARDS = [f"{i:04d}" for i in range(SCATTER_COUNT)]

# Joint genotyping is only meaningful within a reference, so samples are grouped by the
# reference they were mapped against.
REFERENCE_TO_SAMPLES = {}
for _sample, _reference in SAMPLE_TO_REFERENCE.items():
    REFERENCE_TO_SAMPLES.setdefault(_reference, []).append(_sample)

####
# Workflow rules
####

# Several output paths nest a sample or reference directory above a shard file, so the
# wildcards are constrained to stop them from matching across path separators.
wildcard_constraints:
    sample="[^/]+",
    reference="[^/]+",
    shard="[0-9]{4}"

# Build the GATK sequence dictionary for the host reference. The .fna.fai companion is
# produced by rule reference_faidx in preparing.smk, which this branch relies on rather than
# redefining; htsjdk locates the FASTA index by appending .fai to the full filename and the
# dictionary by replacing the extension, hence the two naming conventions.

rule genomics_reference_dict:
    input:
        reference=f"{OUTPUT_DIR}/data/references/{{reference}}.fna",
        fai=f"{OUTPUT_DIR}/data/references/{{reference}}.fna.fai"
    output:
        dict = f"{OUTPUT_DIR}/data/references/{{reference}}.dict"
    params:
        gatk_module={GATK_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 1024 * 10)) * 2 ** (attempt - 1))
    message: "Building sequence dictionary for reference {wildcards.reference}..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module}
        gatk CreateSequenceDictionary -R {input.reference} -O {output.dict}
        """

# The host-mapped alignments arrive from preprocessing already in CRAM and already indexed, so
# this branch goes straight to duplicate marking. That format choice is what makes GATK usable
# here at all: htsjdk reads BAM only through BAI indexes, which cannot address positions beyond
# 2^29-1 (~512 Mbp), whereas the CRAM index stores absolute file offsets and has no such
# ceiling.
#
# Mark (not remove) optical and PCR duplicates. HaplotypeCaller excludes flagged duplicates
# itself, and keeping them in the file preserves the ability to recompute metrics later.

rule genomics_mark_duplicates:
    input:
        cram=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.cram",
        crai=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.cram.crai",
        reference=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.fna",
        fai=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.fna.fai",
        dict=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.dict"
    output:
        cram=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}.cram",
        crai=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}.cram.crai",
        metrics=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}.metrics.txt"
    params:
        samtools_module={SAMTOOLS_MODULE},
        gatk_module={GATK_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE},
        tempdir=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}_tmp"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(16*1024, int(input.size_mb * 3)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(30, int(input.size_mb / 50)) * 2 ** (attempt - 1))
    message: "Marking duplicates in {wildcards.sample}..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module} {params.samtools_module}
        mkdir -p {params.tempdir}
        gatk --java-options "-Xmx$(({resources.mem_mb} * 8 / 10))m" MarkDuplicates \
            -I {input.cram} \
            -O {output.cram} \
            -M {output.metrics} \
            -R {input.reference} \
            --TMP_DIR {params.tempdir} \
            --CREATE_INDEX false
        samtools index -@ {threads} {output.cram}
        rm -rf {params.tempdir}
        """

# Split the reference into evenly sized interval shards so variant calling can be parallelised.
# The default subdivision mode splits large contigs, which is what makes a multi-Gbp chromosome
# tractable — BALANCING_WITHOUT_INTERVAL_SUBDIVISION would keep it in a single shard.

rule genomics_split_intervals:
    input:
        reference=f"{OUTPUT_DIR}/data/references/{{reference}}.fna",
        fai=f"{OUTPUT_DIR}/data/references/{{reference}}.fna.fai",
        dict=f"{OUTPUT_DIR}/data/references/{{reference}}.dict"
    output:
        intervals=expand(
            f"{OUTPUT_DIR}/genomics/intervals/{{{{reference}}}}/{{shard}}-scattered.interval_list",
            shard=SHARDS
        )
    params:
        gatk_module={GATK_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE},
        scatter_count=SCATTER_COUNT,
        outdir=f"{OUTPUT_DIR}/genomics/intervals/{{reference}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(15 * 2 ** (attempt - 1))
    message: "Splitting reference {wildcards.reference} into {params.scatter_count} interval shards..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module}
        gatk SplitIntervals \
            -R {input.reference} \
            --scatter-count {params.scatter_count} \
            -O {params.outdir}
        """

# Per-sample, per-shard variant calling in GVCF mode.
#
# GATK would normally write a .tbi alongside a bgzipped output, but tabix indexes share the
# 512 Mbp ceiling of BAI and cannot index a multi-Gbp chromosome. Index creation is therefore
# disabled here and a CSI index is built with bcftools instead.

rule genomics_haplotype_caller:
    input:
        cram=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}.cram",
        crai=f"{OUTPUT_DIR}/genomics/markdup/{{sample}}.cram.crai",
        reference=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.fna",
        fai=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.fna.fai",
        dict=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}.dict",
        intervals=lambda wildcards: f"{OUTPUT_DIR}/genomics/intervals/{SAMPLE_TO_REFERENCE[wildcards.sample]}/{wildcards.shard}-scattered.interval_list"
    output:
        gvcf=f"{OUTPUT_DIR}/genomics/gvcf/{{sample}}/{{shard}}.g.vcf.gz",
        csi=f"{OUTPUT_DIR}/genomics/gvcf/{{sample}}/{{shard}}.g.vcf.gz.csi"
    params:
        gatk_module={GATK_MODULE},
        bcftools_module={BCFTOOLS_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE}
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(16*1024, int(input.size_mb * 2)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(60, int(input.size_mb / 10)) * 2 ** (attempt - 1))
    message: "Calling variants in {wildcards.sample} (shard {wildcards.shard})..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module} {params.bcftools_module}
        gatk --java-options "-Xmx$(({resources.mem_mb} * 8 / 10))m" HaplotypeCaller \
            -R {input.reference} \
            -I {input.cram} \
            -L {input.intervals} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            --create-output-variant-index false
        bcftools index -c {output.gvcf}
        """

# Combine the per-sample GVCFs of a shard, then joint-genotype them.
#
# CombineGVCFs scales linearly with cohort size; if this branch ever runs on a large cohort,
# GenomicsDBImport is the drop-in replacement here.

rule genomics_combine_gvcfs:
    input:
        reference=f"{OUTPUT_DIR}/data/references/{{reference}}.fna",
        fai=f"{OUTPUT_DIR}/data/references/{{reference}}.fna.fai",
        dict=f"{OUTPUT_DIR}/data/references/{{reference}}.dict",
        intervals=f"{OUTPUT_DIR}/genomics/intervals/{{reference}}/{{shard}}-scattered.interval_list",
        gvcfs=lambda wildcards: expand(
            f"{OUTPUT_DIR}/genomics/gvcf/{{sample}}/{wildcards.shard}.g.vcf.gz",
            sample=REFERENCE_TO_SAMPLES[wildcards.reference]
        ),
        csis=lambda wildcards: expand(
            f"{OUTPUT_DIR}/genomics/gvcf/{{sample}}/{wildcards.shard}.g.vcf.gz.csi",
            sample=REFERENCE_TO_SAMPLES[wildcards.reference]
        )
    output:
        gvcf=temp(f"{OUTPUT_DIR}/genomics/combined/{{reference}}/{{shard}}.g.vcf.gz"),
        csi=temp(f"{OUTPUT_DIR}/genomics/combined/{{reference}}/{{shard}}.g.vcf.gz.csi")
    params:
        gatk_module={GATK_MODULE},
        bcftools_module={BCFTOOLS_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE},
        variant_args=lambda wildcards, input: " ".join(f"-V {gvcf}" for gvcf in input.gvcfs)
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(16*1024, int(input.size_mb * 3)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(60, int(input.size_mb / 10)) * 2 ** (attempt - 1))
    message: "Combining GVCFs of reference {wildcards.reference} (shard {wildcards.shard})..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module} {params.bcftools_module}
        gatk --java-options "-Xmx$(({resources.mem_mb} * 8 / 10))m" CombineGVCFs \
            -R {input.reference} \
            -L {input.intervals} \
            {params.variant_args} \
            -O {output.gvcf} \
            --create-output-variant-index false
        bcftools index -c {output.gvcf}
        """

rule genomics_genotype_gvcfs:
    input:
        reference=f"{OUTPUT_DIR}/data/references/{{reference}}.fna",
        fai=f"{OUTPUT_DIR}/data/references/{{reference}}.fna.fai",
        dict=f"{OUTPUT_DIR}/data/references/{{reference}}.dict",
        intervals=f"{OUTPUT_DIR}/genomics/intervals/{{reference}}/{{shard}}-scattered.interval_list",
        gvcf=f"{OUTPUT_DIR}/genomics/combined/{{reference}}/{{shard}}.g.vcf.gz",
        csi=f"{OUTPUT_DIR}/genomics/combined/{{reference}}/{{shard}}.g.vcf.gz.csi"
    output:
        vcf=temp(f"{OUTPUT_DIR}/genomics/genotyped/{{reference}}/{{shard}}.vcf.gz"),
        csi=temp(f"{OUTPUT_DIR}/genomics/genotyped/{{reference}}/{{shard}}.vcf.gz.csi")
    params:
        gatk_module={GATK_MODULE},
        bcftools_module={BCFTOOLS_MODULE},
        java_module={JAVA_MODULE},
        python_module={PYTHON_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(16*1024, int(input.size_mb * 3)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(60, int(input.size_mb / 10)) * 2 ** (attempt - 1))
    message: "Genotyping reference {wildcards.reference} (shard {wildcards.shard})..."
    shell:
        """
        module purge
        module load {params.python_module} {params.java_module} {params.gatk_module} {params.bcftools_module}
        gatk --java-options "-Xmx$(({resources.mem_mb} * 8 / 10))m" GenotypeGVCFs \
            -R {input.reference} \
            -L {input.intervals} \
            -V {input.gvcf} \
            -O {output.vcf} \
            --create-output-variant-index false
        bcftools index -c {output.vcf}
        """

# Concatenate the shards back into a single cohort VCF. bcftools concat is used rather than
# GatherVcfs because it writes a CSI index directly and never touches a tabix index.

rule genomics_gather_vcfs:
    input:
        vcfs=expand(
            f"{OUTPUT_DIR}/genomics/genotyped/{{{{reference}}}}/{{shard}}.vcf.gz",
            shard=SHARDS
        ),
        csis=expand(
            f"{OUTPUT_DIR}/genomics/genotyped/{{{{reference}}}}/{{shard}}.vcf.gz.csi",
            shard=SHARDS
        )
    output:
        vcf=f"{OUTPUT_DIR}/genomics/variants/{{reference}}.vcf.gz",
        csi=f"{OUTPUT_DIR}/genomics/variants/{{reference}}.vcf.gz.csi"
    params:
        bcftools_module={BCFTOOLS_MODULE}
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 2)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(30, int(input.size_mb / 50)) * 2 ** (attempt - 1))
    message: "Gathering cohort variants for reference {wildcards.reference}..."
    shell:
        """
        module purge
        module load {params.bcftools_module}
        bcftools concat --allow-overlaps --threads {threads} -Oz -o {output.vcf} {input.vcfs}
        bcftools index -c {output.vcf}
        """
