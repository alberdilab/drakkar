####
# Define config variables
####

DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

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
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against genome catalogue..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

rule quantify_reads:
    input:
        expand(f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam", sample=samples)
    output:
        count_table = "3_Outputs/10_Final_tables/unfiltered_count_table.txt",
        mapping_rate = "3_Outputs/10_Final_tables/mapping_rate.txt"
    params:
        BAMs = "3_Outputs/9_MAG_catalogue_mapping/BAMs",
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=64,
        time='02:00:00'
    benchmark:
        "3_Outputs/0_Logs/coverm.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/coverm.log"
    message:
        "Creating the count table with CoverM"
    shell:
        """
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_fraction length \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.count_table}

        #relative abundance for report
        coverm genome \
            -b {params.BAMs}/*.bam \
            -s ^ \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.mapping_rate}
        """

rule gtdb:

rule kegg:

rule gems:
