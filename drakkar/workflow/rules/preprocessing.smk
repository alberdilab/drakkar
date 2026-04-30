####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
FASTP_MODULE = config["FASTP_MODULE"]
SINGLEM_MODULE = config["SINGLEM_MODULE"]

####
# Workflow rules
####

# Reference preparation, mapping, metagenomic read extraction, and host read summaries only run
# if a reference genome or prebuilt reference index is provided.

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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
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

# Run SingleM on the filtered read pairs to generate OTU tables and condensed
# taxonomic profiles for microbial fraction estimation.

rule singlem:
    input:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
    output:
        otu=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_OTU.tsv",
        condense=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_cond.tsv"
    params:
        singlem_module={SINGLEM_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.singlem_module}
        singlem pipe \
            -1 {input.r1} \
            -2 {input.r2} \
            --otu-table {output.otu} \
            --taxonomic-profile {output.condense} \
            --threads {threads}
        """

# Estimate microbial fraction for each sample using SingleM.

rule singlem_mf:
    input:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz",
        profile=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_cond.tsv"
    output:
        f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_smf.tsv"
    params:
        singlem_module={SINGLEM_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.singlem_module}
        singlem microbial_fraction \
            -1 {input.r1} \
            -2 {input.r2} \
            --input-profile {input.profile} \
            --output-tsv {output}
        """

rule preprocessings_stats:
    input:
        fastp=expand(f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json", sample=samples),
        fastq=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz", sample=samples),
        singlem=expand(f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_smf.tsv", sample=samples) if FRACTION else [],
        nonpareil=expand(f"{OUTPUT_DIR}/preprocessing/nonpareil/{{sample}}_np.tsv", sample=samples) if NONPAREIL else []
    output:
        f"{OUTPUT_DIR}/preprocessing.tsv"
    localrule: True
    params:
        package_dir={PACKAGE_DIR},
        singlem_arg=lambda wildcards, input: "-s " + " ".join(input.singlem) if input.singlem else "",
        nonpareil_arg=lambda wildcards, input: "-n " + " ".join(input.nonpareil) if input.nonpareil else ""
    threads: 1
    resources:
        mem_mb= 1*1024,
        runtime= 5
    message: "Creating preprocessing stats..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/preprocessing_stats.py -p {input.fastp} -f {input.fastq} {params.singlem_arg} {params.nonpareil_arg} -o {output}
        """
