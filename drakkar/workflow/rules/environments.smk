# Toy pipeline that only declares environments
PACKAGE_DIR = config["package_dir"]

rule all:
    input:
        f"drakkar_environments/env4.txt"

rule env1:
    output:
        f"drakkar_environments/env1.txt"
    threads: 1
    conda: 
        f"{PACKAGE_DIR}/workflow/envs/cataloging.yaml"
    shell:
        """
        touch drakkar_environments/env1.txt
        """
    
rule env2:
    input:
        f"drakkar_environments/env1.txt"
    output:
        f"drakkar_environments/env2.txt"
    threads: 1
    conda: 
        f"{PACKAGE_DIR}/workflow/envs/profiling_genomes.yaml"
    shell:
        """
        touch drakkar_environments/env2.txt
        """

rule env3:
    input:
        f"drakkar_environments/env2.txt"
    output:
        f"drakkar_environments/env3.txt"
    threads: 1
    conda: 
        f"{PACKAGE_DIR}/workflow/envs/annotating_taxonomu.yaml"
    shell:
        """
        touch drakkar_environments/env3.txt
        """

rule env4:
    input:
        f"drakkar_environments/env3.txt"
    output:
        f"drakkar_environments/env4.txt"
    threads: 1
    conda: 
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    shell:
        """
        touch drakkar_environments/env4.txt
        """