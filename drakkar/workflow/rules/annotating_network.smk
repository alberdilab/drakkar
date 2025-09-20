####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
EMAPPER_MODULE = config["EMAPPER_MODULE"]
PATHWAYTOOLS_MODULE = config["PATHWAYTOOLS_MODULE"]
BLAST_MODULE = config["BLAST_MODULE"]
M2M_MODULE = config["M2M_MODULE"]
GAPSEQ_MODULE = config["GAPSEQ_MODULE"]

# Annotation databases
EGGNOG_DB = config["EGGNOG_DB"]
CARVEME_DB = config["CARVEME_DB"]
GAPSEQ_DB = config["GAPSEQ_DB"]

####
# Workflow rules
####

rule prodigal:
    input:
        lambda wildcards: MAGS_TO_FILES[wildcards.mag]
    output:
        fna=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.fna",
        faa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa",
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff"
    log:
        f"{OUTPUT_DIR}/log/annotating/prodigal/{{mag}}.log"
    params:
         prodigal_module={PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    threads: 1
    shell:
        """
        module load {params.prodigal_module}
        set -euo pipefail
        mkdir -p "$(dirname {output.gff})"
        if [[ "{input}" == *.gz ]]; then
            echo "Running gzip"  
            gzip -dc {input} | prodigal -i - -d {output.fna} -a {output.faa} -o {output.gff} -f gff
        else
            echo "Running non-gzip"  
            prodigal -i {input} -d {output.fna} -a {output.faa} -o {output.gff} -f gff
        fi
        """

rule carveme:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        f"{OUTPUT_DIR}/annotating/carveme/{{mag}}.sbml"
    params:
        carveme_db={CARVEME_DB}
    threads:
        8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 8) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 100) * 2 ** (attempt - 1))
    shell:
        """
        carve {input} -o {output} --mediadb {params.carveme_db} --gapfill M3 --fbc2
        """

rule gapseq:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}.xml"
    params:
        gapseq_db={GAPSEQ_DB},
        gapseq_dir=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}"
    threads:
        8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 8) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 500) * 2 ** (attempt - 1))
    shell:
        """
        set -euo pipefail
        mkdir -p {params.gapseq_dir}
        cd {params.gapseq_dir}
        gapseq doall {input} {params.gapseq_db}/gut.csv -K {threads}
        """

rule mergem:
    input:
        carveme=f"{OUTPUT_DIR}/annotating/carveme/{{mag}}.sbml",
        gapseq=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}.xml"
    output:
        sbml=f"{OUTPUT_DIR}/annotating/mergem/{{mag}}.sbml",
        report=f"{OUTPUT_DIR}/annotating/mergem/{{mag}}.tsv"
    threads:
        1
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 8) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 300) * 2 ** (attempt - 1))
    shell:
        """
        mergem -i {input.carveme} {input.gapseq} -o {output.sbml} --report {output.report}
        """
