####
# Define config variables
####

FASTP_MODULE = config["FASTP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

####
# Calculate file sizes
####

reads_mb = calculate_file_sizes(READS_DIR)
reads_mb = {key.replace('_1.fq.gz', ''): value for key, value in reads_mb.items()}
reads_mb_total = sum(reads_mb.values())

####
# Run preprocessing rules
####

# The rules reference_index, reference_map, metagenomic_reads and host_reads are only run if
# the reference genome file is provided. Otherwise, the fastp rule already outputs the final files.

if REFERENCE:

    reference_mb = calculate_file_size(REFERENCE)

    rule fastp:
        input:
            r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
        output:
            r1=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_2.fq.gz",
            html=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.html",
            json=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json"
        params:
            fastp_module={FASTP_MODULE}
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 5) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.fastp_module}
            fastp \
                --in1 {input.r1} --in2 {input.r2} \
                --out1 {output.r1} --out2 {output.r2} \
                --trim_poly_g \
                --trim_poly_x \
                --low_complexity_filter \
                --n_base_limit 5 \
                --qualified_quality_phred 20 \
                --length_required 60 \
                --thread {threads} \
                --html {output.html} \
                --json {output.json} \
                --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
            """

    rule reference_index:
        input:
            f"{OUTPUT_DIR}/data/references/{{reference}}.fna"
        output:
            index=f"{OUTPUT_DIR}/data/references/{{reference}}.rev.1.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=lambda wildcards: f"{OUTPUT_DIR}/data/references/{wildcards.reference}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(1,int(reference_mb * 10 * 2 ** (attempt - 1))),
            runtime=lambda wildcards, attempt: max(15,int((reference_mb / 20 * 2 ** (attempt - 1))))
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            cat {input} > {params.basename}.fna
            """

    rule reference_map:
        input:
            index=lambda wildcards: expand(
                f"{OUTPUT_DIR}/data/references/{{reference}}.rev.1.bt2",
                reference=[SAMPLE_TO_REFERENCE[wildcards.sample]]
            ),
            r1=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=lambda wildcards: f"{OUTPUT_DIR}/data/references/{SAMPLE_TO_REFERENCE[wildcards.sample]}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * reference_mb / 500) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * reference_mb / 10) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
            """

    rule metagenomic_reads:
        input:
            f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
        output:
            r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE}
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            samtools view -b -f12 -@ {threads} {input} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -
            """

    rule host_reads:
        input:
            f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
        output:
            f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.bam"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE}
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 10) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            samtools view -b -F12 -@ {threads} {input} | samtools sort -@ {threads} -o {output} -
            """
else:
    rule fastp:
        input:
            r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
        output:
            r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
            r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz",
            html=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.html",
            json=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json"
        params:
            fastp_module={FASTP_MODULE}
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 5) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.fastp_module}
            fastp \
                --in1 {input.r1} --in2 {input.r2} \
                --out1 {output.r1} --out2 {output.r2} \
                --trim_poly_g \
                --trim_poly_x \
                --low_complexity_filter \
                --n_base_limit 5 \
                --qualified_quality_phred 20 \
                --length_required 60 \
                --thread {threads} \
                --html {output.html} \
                --json {output.json} \
                --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
            """
