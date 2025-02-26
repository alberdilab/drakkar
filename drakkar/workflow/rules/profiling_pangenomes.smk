####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]

####
# Workflow rules
####

# Predict genes from all MAGs
rule prodigal:
    input:
        lambda wildcards: BINS_TO_FILES[wildcards.bin]
    output:
        nt=f"{OUTPUT_DIR}/profiling_pangenomes/prodigal/{{bin}}.fna",
        aa=f"{OUTPUT_DIR}/profiling_pangenomes/prodigal/{{bin}}.faa"
    params:
        prodigal_module={PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 3) * 2 ** (attempt - 1))
    threads: 1
    shell:
        """
        module load {params.prodigal_module}
        prodigal -i {input} -d {output.nt} -a {output.aa}
        """

# Create a reference table with genome and gene names to trace origin of genes
rule create_genome_gene_table:
    input:
        expand(f"{OUTPUT_DIR}/profiling_pangenomes/prodigal/{{bin}}.fna", bin=bins)
    output:
        f"{OUTPUT_DIR}/profiling_pangenomes/genome_gene.csv"
    shell:
        """
        python {params.package_dir}/workflow/scripts/genome_gene_table.py {input} {output}
        """

checkpoint dereplicate:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values()),
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    output:
        Cdb=f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv",
        Bdb=f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Bdb.csv",
        Wdb=f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Wdb.csv"
    params:
        drep_module={DREP_MODULE},
        outdir=f"{OUTPUT_DIR}/profiling_pangenomes/drep/"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 20) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.drep_module}
        rm -rf {params.outdir}
        dRep compare {params.outdir} -p {threads} -g {input.genomes} -sa 0.95 --genomeInfo {input.metadata}
        """

# Cluster original genomes into group of genomes to be dereplicated using mmseqs
# Note this this a checkpoint. The outputs of this rule cannot be predicted, as
# they depend on the content of the Cdb.csv file generated in the previous drep rule.
# In consequence, the checkpoint forces snakemake to re-create the pipeline structure
# and consider the new wildards. Note the use of functions to generate the dynamic
# inputs for the direct dependencies of this rule in the subsequent steps.

checkpoint cluster_genomes:
    input:
        csv=f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv",
        nt=expand(f"{OUTPUT_DIR}/profiling_pangenomes/prodigal/{{bin}}.fna", bin=bins),
        aa=expand(f"{OUTPUT_DIR}/profiling_pangenomes/prodigal/{{bin}}.faa", bin=bins)
    output:
        directory(f"{OUTPUT_DIR}/profiling_pangenomes/clusters")
    run:
        import pandas as pd
        import os

        # Read clustering information from the CSV
        df = pd.read_csv(input.csv)
        cluster_dict = df.groupby("secondary_cluster")["genome"].apply(list).to_dict()

        os.makedirs(CLUSTER_DIR, exist_ok=True)

        # Create cluster files by concatenating genomes
        for cluster_id, genome_list in cluster_dict.items():
            # Generate nucleotide sequences
            cluster_file_nt = f"{OUTPUT_DIR}/profiling_pangenomes/clusters/cluster_{cluster_id}.fna"
            with open(cluster_file_nt, "w") as out_file:
                for genome_id in genome_list:
                    prodigal_file = f"{PRODIGAL_DIR}/{genome_id}"
                    if prodigal_file in input.nt:
                        with open(prodigal_file, "r") as f:
                            out_file.write(f.read())
            # Generate amino acid sequences
            cluster_file_aa = f"{OUTPUT_DIR}/profiling_pangenomes/clusters/cluster_{cluster_id}.faa"
            with open(cluster_file_aa, "w") as out_file:
                for genome_id in genome_list:
                    genome_id_aa = genome_id.replace(".fna", ".faa")
                    prodigal_file_aa = f"{PRODIGAL_DIR}/{genome_id_aa}"
                    if prodigal_file_aa in input.aa:
                        with open(prodigal_file_aa, "r") as f:
                            out_file.write(f.read())

# Functions to define the input files dynamically.

def get_cluster_ids_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    return df["secondary_cluster"].unique()

## These functions return each file separately, and need to be combined with a lambda expression in the input

def get_cluster_fna_sep(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return f"{CLUSTER_DIR}/cluster_{wildcards.cluster_id}.fna"

def get_gene_fna_sep(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return f"{MMSEQS_DIR}/cluster_{wildcards.cluster_id}.fna"

def get_gene_faa_sep(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return f"{MMSEQS_DIR}/cluster_{wildcards.cluster_id}.faa"

def get_gene_tsv_sep(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return f"{MMSEQS_DIR}/cluster_{wildcards.cluster_id}.tsv"

## These functions return a list of all desired files, and can be called directly in the input

def get_gene_fna(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return expand(f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}.fna", cluster_id=cluster_ids)

def get_gene_faa(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return expand(f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}.faa", cluster_id=cluster_ids)

def get_gene_tsv(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return expand(f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}.tsv", cluster_id=cluster_ids)

def get_kofams_tsv(wildcards):
    checkpoint_output = checkpoints.cluster_genomes.get(**wildcards).output[0]
    cluster_ids = get_cluster_ids_from_csv(f"{OUTPUT_DIR}/profiling_pangenomes/drep/data_tables/Cdb.csv")
    return expand(f"{KOFAMS_DIR}/cluster_{{cluster_id}}.tsv", cluster_id=cluster_ids)

####
#### The rules below are only executed after the checkpoint
####

# Dereplicate the genes to obtain a list of reference genes for the cluster
rule run_mmseqs:
    input:
        cluster_genes=lambda wildcards: get_cluster_fna_sep(wildcards)
    output:
        fna=f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}.fna",
        tsv=f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}.tsv"
    params:
        tmp=f"{OUTPUT_DIR}/profiling_pangenomes/mmseqs/cluster_{{cluster_id}}_tmp"
    conda:
        "environments/clustering.yml"
    shell:
        """
        mmseqs easy-linclust {input.cluster_genes} {MMSEQS_DIR}/cluster_{wildcards.cluster_id}_output {params.tmp} --min-seq-id 0.8 --cov-mode 1 -c 0.8
        rm {MMSEQS_DIR}/cluster_{wildcards.cluster_id}_output_all_seqs.fasta
        mv {MMSEQS_DIR}/cluster_{wildcards.cluster_id}_output_rep_seq.fasta {output.fna}
        mv {MMSEQS_DIR}/cluster_{wildcards.cluster_id}_output_cluster.tsv {output.tsv}
        rm -rf {params.tmp}
        """
