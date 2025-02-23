####
# Define config variables
####

DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

rule dereplicate:
    input:
         expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        dir=directory(f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes"),
        Cdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Cdb.csv",
        Bdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Bdb.csv",
        Wdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Wdb.csv"
    params:
        drep_module={DREP_MODULE}
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 10) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        source deactivate
        module load {params.drep_module}
        dRep dereplicate {output.dir} -p {threads} -g {input}
        """

rule merge_catalogue:
    input:
        directory(f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes")
    output:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    message: "Merging genomes into a single catalogue..."
    shell:
        """
        cat {input}/* > {output}
        """

rule index_catalogue:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    output:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        """

rule map_to_catalogue:
    input:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2",
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
    output:
        f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(reads_mb.get(wildcards.sample, 1) * 2 * 1024 / 500) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(reads_mb.get(wildcards.sample, 1) / 1024 *  1024 / 10) * 2 ** (attempt - 1))
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
