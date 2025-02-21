####
# Define config variables
####

PACKAGE_DIR = config["package_dir"]
FASTP_MODULE = config["FASTP_MODULE"]

####
# Run preprocessing rules
####

# The rules reference_index, reference_map, metagenomic_reads and host_reads are only run if
# the reference genome file is provided. Otherwise, the fastp rule already outputs the final files.

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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024) * 2 ** (attempt - 1))
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

rule preprocessings_stats:
    input:
        fastp=expand(f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json", sample=samples),
        fastq=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz", sample=samples),
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
        python {params.package_dir}/workflow/scripts/preprocessing_stats.py -p "{input.fastp}" -f "{input.fastq}" -o {output}
        """
