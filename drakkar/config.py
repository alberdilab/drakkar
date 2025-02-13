import os

# Define software modules
SNAKEMAKE_MODULE = os.getenv("SNAKEMAKE_MODULE", "snakemake/8.16.0")
FASTP_MODULE = os.getenv("FASTP_MODULE", "fastp/0.23.4")
