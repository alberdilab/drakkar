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

rule gapseq_find:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        pathways=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-all-Pathways.tbl",
        reactions=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-all-Reactions.tbl"
    params:
        gapseq_db={GAPSEQ_DB},
        gapseq_dir=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}"
    threads:
        8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 500) * 2 ** (attempt - 1))
    shell:
        """
        set -euo pipefail
        mkdir -p {params.gapseq_dir}
        cd {params.gapseq_dir}
        gapseq find -p all -K {threads} {input} 
        """

rule gapseq_find_transport:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        transporter=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-Transporter.tbl"
    params:
        gapseq_db={GAPSEQ_DB},
        gapseq_dir=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}"
    threads:
        8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 25) * 2 ** (attempt - 1))
    shell:
        """
        set -euo pipefail
        mkdir -p {params.gapseq_dir}
        cd {params.gapseq_dir}
        gapseq find-transport -K {threads} {input} 
        """

rule gapseq_draft:
    input:
        fasta=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa",
        pathways=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-all-Pathways.tbl",
        reactions=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-all-Reactions.tbl",
        transporter=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-Transporter.tbl"
    output:
        model=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-draft.RDS",
        weights=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-rxnWeights.RDS",
        Xgenes=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-rxnXgenes.RDS"
    params:
        gapseq_db={GAPSEQ_DB},
        gapseq_dir=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}"
    threads:
        8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/annotating_network.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 600) * 2 ** (attempt - 1))
    shell:
        """
        set -euo pipefail
        mkdir -p {params.gapseq_dir}
        cd {params.gapseq_dir}
        gapseq draft -r {input.reactions} -t {input.transporter} -p {input.pathways} -c {input.fasta} -K {threads}
        """

rule gapseq_fill:
    input:
        model=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-draft.RDS",
        weights=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-rxnWeights.RDS",
        Xgenes=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}/{{mag}}-rxnXgenes.RDS"
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 500) * 2 ** (attempt - 1))
    shell:
        """
        set -euo pipefail
        cd {params.gapseq_dir}
        gapseq fill -m {input.model} -c {input.weights} -g {input.Xgenes} -n {params.gapseq_db}/gut.csv
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
