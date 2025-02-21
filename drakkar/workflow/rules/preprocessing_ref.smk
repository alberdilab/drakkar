####
# Define config variables
####

PACKAGE_DIR = config["package_dir"]
FASTP_MODULE = config["FASTP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]


####
# Run preprocessing rules
####

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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 3) * 2 ** (attempt - 1))
    message: "Quality-filtering sample {wildcards.sample}..."
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
        basename=f"{OUTPUT_DIR}/data/references/{{reference}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Indexing reference genome {wildcards.reference}..."
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
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against reference genome..."
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
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz",
        reads=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metareads",
        bases=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metabases"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Extracting metagenomic reads of {wildcards.sample}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        samtools view -b -f12 -@ {threads} {input} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -
        samtools view -b -f12 -@ {threads} {input} | samtools view -c {input} > {output.reads}
        samtools view -b -f12 -@ {threads} {input} | samtools stats - | grep "^SN" | grep "bases mapped (cigar)" | cut -f3 > {output.bases}
        """

rule host_reads:
    input:
        f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
    output:
        bam=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.bam",
        reads=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostreads",
        bases=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostbases"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Extracting host reads of {wildcards.sample}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        samtools view -b -F12 -@ {threads} {input} | samtools sort -@ {threads} -o {output.bam} -
        samtools view -b -F12 -@ {threads} {input} | samtools view -c {input} > {output.reads}
        samtools view -b -F12 -@ {threads} {input} | samtools stats - | grep "^SN" | grep "bases mapped (cigar)" | cut -f3 > {output.bases}
        """

rule preprocessings_stats:
    input:
        fastp=expand(f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json", sample=samples),
        metagenomic=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz", sample=samples),
        genomic=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.bam", sample=samples)
    output:
        f"{OUTPUT_DIR}/preprocessing.tsv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    message: "Creating preprocessing stats..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/preprocessing_stats.py -f "{input.fastp}" -m "{input.metagenomic}" -g "{input.genomic}" -o {output}
        """
