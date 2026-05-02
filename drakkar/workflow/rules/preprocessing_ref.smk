####
# Variables parsed from the config.yaml file
####

PACKAGE_DIR = config["package_dir"]

# Software modules
FASTP_MODULE = config["FASTP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
SINGLEM_MODULE = config["SINGLEM_MODULE"]
SINGLEM_DB = config["SINGLEM_DB"]
MULTIQC_MODULE = config["MULTIQC_MODULE"]

####
# Workflow rules
####

# Run fastp on paired-end reads to perform adapter/quality trimming and filtering.  
# Outputs the cleaned FASTQ files along with an HTML and JSON report that MultiQC can parse.

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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 1024 * 3) * 2 ** (attempt - 1)))
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

# Align quality-filtered paired-end reads to the corresponding reference genome using Bowtie2.  
# The alignments are converted to BAM format and sorted with samtools for downstream analyses.

rule reference_map:
    input:
        index=lambda wildcards: expand(
            f"{OUTPUT_DIR}/data/references/{{reference}}.index.ready",
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 20) * 2 ** (attempt - 1)))
    message: "Mapping {wildcards.sample} against reference genome..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

# Generate alignment quality metrics for each BAM file using samtools.  
# Produces an index (.bai), a flagstat summary, idxstats, and detailed stats file for MultiQC reporting.

rule samtools_stats:
    input:
        rules.reference_map.output
    output:
        bai      = f"{OUTPUT_DIR}/preprocessing/samtools/{{sample}}.bam.bai",
        flagstat = f"{OUTPUT_DIR}/preprocessing/samtools/{{sample}}.flagstat.txt",
        idxstats = f"{OUTPUT_DIR}/preprocessing/samtools/{{sample}}.idxstats.txt",
        stats    = f"{OUTPUT_DIR}/preprocessing/samtools/{{sample}}.stats.txt"
    params:
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 2) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(10 * 2 ** (attempt - 1))
    message: "Generating mapping stats for {wildcards.sample}..."
    shell:
        """
        module load {params.samtools_module}
        samtools index {input} {output.bai}
        samtools flagstat {input} > {output.flagstat}
        samtools idxstats {input} > {output.idxstats}
        samtools stats {input} > {output.stats}
        """

# Split mapped BAM files into metagenomic (unmapped) and host (mapped) read sets.
# Count sidecars are temporary inputs for preprocessing.tsv and are removed after use.

rule split_reads:
    input:
        f"{OUTPUT_DIR}/preprocessing/bowtie2/{{sample}}.bam"
    output:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz",
        metareads=temp(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metareads"),
        metabases=temp(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metabases"),
        bam=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.bam",
        hostreads=temp(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostreads"),
        hostbases=temp(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostbases")
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Extracting metagenomic reads of {wildcards.sample}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        samtools view -b -f12 -@ {threads} {input} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -
        samtools view -b -f12 -@ {threads} {input} | samtools view -c - > {output.metareads}
        samtools view -f12 -@ {threads} {input} | awk '{{sum += length($10)}} END {{print sum}}' > {output.metabases}
        samtools view -b -F12 -@ {threads} {input} | samtools sort -@ {threads} -o {output.bam} -
        samtools view -b -F12 -@ {threads} {input} | samtools view -c - > {output.hostreads}
        samtools view -F12 -@ {threads} {input} | awk '{{sum += length($10)}} END {{print sum}}' > {output.hostbases}
        """

# Run SingleM on the filtered metagenomic read pairs to generate OTU tables  
# and condensed taxonomic profiles for downstream community analysis.

rule singlem:
    input:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
    output:
        otu=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_OTU.tsv",
        condense=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_cond.tsv"
    params:
        singlem_module={SINGLEM_MODULE},
        singlem_db={SINGLEM_DB}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    shell:
        """
        module load {params.singlem_module}
        singlem pipe \
            -1 {input.r1} \
            -2 {input.r2} \
            --otu-table {output.otu} \
            --taxonomic-profile {output.condense} \
            --metapackage {params.singlem_db} \
            --threads {threads}
        """

# Estimate microbial fraction for each sample using SingleM.  
# Combines paired reads with the SingleM taxonomic profile to produce a microbial fraction table.

rule singlem_mf:
    input:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz",
        profile=f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_cond.tsv"
    output:
        f"{OUTPUT_DIR}/preprocessing/singlem/{{sample}}_smf.tsv"
    params:
        singlem_module={SINGLEM_MODULE},
        singlem_db={SINGLEM_DB}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    shell:
        """
        module load {params.singlem_module}
        singlem microbial_fraction \
            -1 {input.r1} \
            -2 {input.r2} \
            --input-profile {input.profile} \
            --metapackage {params.singlem_db} \
            --output-tsv {output}
        """

# Aggregate preprocessing statistics across all samples.  
# Combines fastp JSON reports with host/metagenomic read and base counts  
# to produce a single summary table (preprocessing.tsv).

rule preprocessing_stats:
    input:
        fastp=expand(f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json", sample=samples),
        reads_metagenomic=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metareads", sample=samples),
        bases_metagenomic=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.metabases", sample=samples),
        reads_host=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostreads", sample=samples),
        bases_host=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}.hostbases", sample=samples),
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
        mem_mb=cap_mem_mb(1*1024),
        runtime=cap_runtime(5)
    message: "Creating preprocessing stats..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/preprocessing_stats.py -p {input.fastp} -m {input.bases_metagenomic} -M {input.reads_metagenomic} -g {input.bases_host} -G {input.reads_host} {params.singlem_arg} {params.nonpareil_arg} -o {output}
        """

# Generate an HTML report summarizing preprocessing statistics.  
# Uses the aggregated preprocessing table (preprocessing.tsv) and writes a completion flag.
# PROBABLY TO BE DEPRECATED AFTER MULTIQC IMPLEMENTATION

rule preprocessing_report:
    input:
        data=f"{OUTPUT_DIR}/preprocessing.tsv",
        report=f"{OUTPUT_DIR}/drakkar_report.html"
    output:
        done=f"{OUTPUT_DIR}/data/preprocessing_report.done"
    localrule: True
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=cap_mem_mb(1*1024),
        runtime=cap_runtime(5)
    message: "Creating preprocessing stats..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/preprocessing_report.py -i {input.data} -r {input.report} -o {input.report}
        touch {output.done}
        """

# Run MultiQC on preprocessing outputs (fastp JSON and samtools stats).  
# Produces an HTML summary report and a zipped data directory for further inspection.

rule preprocessing_multiqc:
    input:
        fastp=expand(f"{OUTPUT_DIR}/preprocessing/fastp/{{sample}}.json", sample=samples),
        samtools=expand(f"{OUTPUT_DIR}/preprocessing/samtools/{{sample}}.flagstat.txt", sample=samples)
    output:
        html=f"{OUTPUT_DIR}/preprocessing/preprocessing.html"
    params:
        package_dir={PACKAGE_DIR},
        multiqc_module={MULTIQC_MODULE},
        datadir="preprocessing",
        title="Preprocessing report"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Building MultiQC report..."
    shell:
        """
        module load {params.multiqc_module}
        multiqc --zip-data-dir --outdir {params.datadir} --title "{params.title}" --filename {output.html} {params.datadir}
        """
