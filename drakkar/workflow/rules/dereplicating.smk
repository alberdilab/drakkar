####
# Define config variables
####

DREP_MODULE = config["DREP_MODULE"]

rule compare_ani:
    input:
        "{assembly_dir}/assembly.fna"
    output:
        dir=directory(f"{OUTPUT_DIR}/dereplicating/drep/"),
        Cdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Cdb.csv",
        Bdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Bdb.csv",
        Wdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Wdb.csv"
    params:
        drep_module={DREP_MODULE}
    message: "Calculating pairwise average nucleotide identities between bins..."
    shell:
        """
        module load {params.drep_module}
        dRep compare {output.dir} -g {input}
        """

rule dereplicate:
    input:
        fasta="{assembly_dir}/assembly.fna",
        Cdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Cdb.csv",
        Bdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Bdb.csv",
        Wdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Wdb.csv"
    output:
        directory(f"{OUTPUT_DIR}/dereplicating/final/{{drep}}")
    params:
        drep_module={DREP_MODULE},
        threshold=f"{{drep}}"
    message: "Extracting dereplicated bins at {wildcards.drep}%..."
    shell:
        """
        module load {params.drep_module}
        dRep dereplicate {output.dir} --genomeInfo {input.genomeinfo} -sa {params.threshold}
        
        """

rule dereplicate_final:
    input:
        directory(f"{OUTPUT_DIR}/dereplicating/drep/{{drep}}")
    output:
        directory(f"{OUTPUT_DIR}/dereplicating/final/{{drep}}")
    params:
        drep_module={DREP_MODULE},
        threshold=f"{{drep}}"
    message: "Creating final genome catalogue for {wildcards.drep}%..."
    shell:
        """

        """
