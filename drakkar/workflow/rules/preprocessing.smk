FASTP_MODULE = config["FASTP_MODULE"]

rule fastp:
    input:
        r1=f"{READS_DIR}/{{sample}}_1.fq.gz",
        r2=f"{READS_DIR}/{{sample}}_2.fq.gz"
    output:
        r1=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_2.fq.gz",
        html=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.html",
        json=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json"
    params:
        fastp_module={FASTP_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 16) * 2 ** (attempt - 1)),
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
        f"{REFERENCE}"
    output:
        index=expand(f"{OUTPUT_DIR}/preprocessing/reference/reference.{{ext}}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/preprocessing/reference/reference"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(1,int(reference_mb * 5 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: max(15,int((reference_mb / 4 * 2 ** (attempt - 1))))
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        cat {input} > {params.basename}.fna
        """

rule reference_map:
    input:
        index=expand(f"{OUTPUT_DIR}/preprocessing/reference/reference.{{ext}}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
        r1=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}_2.fq.gz"
    output:
        f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/preprocessing/reference/reference"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 16) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """
