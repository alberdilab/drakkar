####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
COVERM_MODULE = config["COVERM_MODULE"]

####
# Workflow rules
####

checkpoint dereplicate:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values()),
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    output:
        Wdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Wdb.csv",
        Cdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Cdb.csv"
    params:
        drep_module={DREP_MODULE},
        outdir=f"{OUTPUT_DIR}/profiling_genomes/drep/"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.drep_module}
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -p {threads} -g {input.genomes} -sa 0.98 --genomeInfo {input.metadata}
        """

# Functions to define the input files dynamically.
def get_mag_ids_from_drep(csv_path):
    df = pd.read_csv(csv_path)
    return df["genome"].unique()

def get_mag_fna(wildcards):
    checkpoint_output = checkpoints.dereplicate.get(**wildcards).output[0]
    selected_bins = get_mag_ids_from_drep(checkpoint_output)
    return expand(f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes/{{bin_id}}", bin_id=selected_bins)

rule merge_catalogue:
    input:
        get_mag_fna
    output:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    localrule: True
    message: "Merging genomes into a single catalogue..."
    shell:
        """
        cat {input} > {output}
        """

rule index_catalogue:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    output:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        """

rule map_to_catalogue:
    input:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2",
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample]
    output:
        f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against genome catalogue..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

rule quantify_reads_catalogue:
    input:
        expand(f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam", sample=samples)
    output:
        f"{OUTPUT_DIR}/profiling_genomes/coverm/coverm.tsv"
    params:
        coverm_module={COVERM_MODULE}
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb / 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 / 5) * 2 ** (attempt - 1))
    message:
        "Generating mapping statistics with..."
    shell:
        """
        module load {params.coverm_module}
        coverm genome \
            -b {input} \
            -s @ \
            -m count covered_bases \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """

rule split_coverm:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/coverm/coverm.tsv"
    output:
        counts=f"{OUTPUT_DIR}/profiling_genomes/final/counts.tsv",
        breadth=f"{OUTPUT_DIR}/profiling_genomes/final/bases.tsv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb = 1*1024,
        runtime = 5
    message: "Generating count and breadth files..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/split_coverm.py {input} {output.counts} {output.breadth}
        """
