####
# Define config variables
####

DREP_MODULE = config["DREP_MODULE"]

rule dereplicate:
    input:
        f"{OUTPUT_DIR}/data/final_bins.txt"
    output:
        dir=directory(f"{OUTPUT_DIR}/dereplicating/drep/"),
        Cdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Cdb.csv",
        Bdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Bdb.csv",
        Wdb=f"{OUTPUT_DIR}/dereplicating/drep/data_tables/Wdb.csv"
    params:
        drep_module={DREP_MODULE}
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.drep_module}
        dRep dereplicate {output.dir} -g {input}
        """
