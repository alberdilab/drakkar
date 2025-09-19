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
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 300) * 2 ** (attempt - 1))
    shell:
        """
        carve {input} -o {output} --mediadb {params.carveme_db} --gapfill M3 --fbc2
        """

rule gapseq:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}.xml"
    params:
        gapseq_module={GAPSEQ_MODULE},
        gapseq_db={GAPSEQ_DB},
        gapseq_dir=f"{OUTPUT_DIR}/annotating/gapseq/{{mag}}"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 8) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 300) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.gapseq_module}
        set -euo pipefail
        mkdir -p {params.gapseq_dir}
        cd {params.gapseq_dir}
        gapseq {input} {params.gapseq_db}/gut.csv -K {threads}
        """

rule emapper:
    input:
        faa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        ann=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}/{{mag}}.emapper.annotations",
        hit=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}/{{mag}}.emapper.hits",
        ort=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}/{{mag}}.emapper.seed_orthologs"
    params:
        emapper_module={EMAPPER_MODULE},
        eggnog_db={EGGNOG_DB},
        outdir=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 8) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 300) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.emapper_module}
        emapper.py  \
            -i {input.faa} \
            --cpu {threads} \
            --data_dir {params.eggnog_db} \
            -o {wildcards.mag} \
            --output_dir {params.outdir} \
            -m diamond --dmnd_ignore_warnings \
            --itype proteins \
            --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 \
            --tax_scope auto --target_orthologs all --go_evidence non-electronic \
            --pfam_realign none
        """

rule emapper2gbk:
    input:
        fna=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.fna",
        faa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa",
        ann=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}/{{mag}}.emapper.annotations"
    output:
        f"{OUTPUT_DIR}/annotating/gbk/{{mag}}/{{mag}}.gbk"
    threads:
        1
    params:
        emapper_module={EMAPPER_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.emapper_module}
        emapper2gbk genes -fn {input.fna} -fp {input.faa} -a {input.ann} -o {output}
        """

rule m2m:
    input:
        expand(f"{OUTPUT_DIR}/annotating/gbk/{{mag}}/{{mag}}.gbk",mag=mags)
    output:
        sbml=f"{OUTPUT_DIR}/annotating/m2m/sbml/{{mag}}.sbml"
    params:
        pathwaytools_module={PATHWAYTOOLS_MODULE},
        blast_module={BLAST_MODULE},
        m2m_module={M2M_MODULE},
        indir=f"{OUTPUT_DIR}/annotating/eggnog/",
        outdir=f"{OUTPUT_DIR}/annotating/m2m/"
    threads:
        24
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 100) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.pathwaytools_module} {params.blast_module} {params.m2m_module}
        m2m recon -g {params.indir} -o {params.outdir} -c {threads}
        """

rule final_network:
    input:
        f"{OUTPUT_DIR}/annotating/m2m/sbml/{{mag}}.sbml"
    output:
        f"{OUTPUT_DIR}/annotating/m2m/{{mag}}.sbml"
    localrule: True
    threads:
        1
    resources:
        mem_gb=1*1024,
        time=1
    shell:
        """
        mv {input} {output}
        """
