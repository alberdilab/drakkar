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
    envmodules:
        {PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    threads: 1
    shell:
        """
        if [[ "{input}" == *.gz ]]; then
            gzip -dc {input} | prodigal -i - -d {output.fna} -a {output.faa} -o {output.gff} -f gff
        else
            prodigal -i {input} -d {output.fna} -a {output.faa} -o {output.gff} -f gff
        fi
        """

rule emapper:
    input:
        faa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        ann=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}.emapper.annotations",
        hit=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}.emapper.hits",
        ort=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}.emapper.seed_orthologs"
    params:
        outdir=f"{OUTPUT_DIR}/annotating/eggnog/",
        tmpdir="tmp"
    threads:
        8
    envmodules:
        {EMAPPER_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 100) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 100) * 2 ** (attempt - 1))
    shell:
        """
        emapper.py  \
            -i {input.faa} \
            --cpu {threads} \
            --data_dir /projects/mjolnir1/data/databases/eggnog-mapper/20230317/ \
            -o {wildcards.mag} \
            --output_dir {params.outdir} \
            --temp_dir {params.tmpdir} \
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
        ann=f"{OUTPUT_DIR}/annotating/eggnog/{{mag}}.emapper.annotations"
    output:
        f"{OUTPUT_DIR}/annotating/gbk/{{mag}}/{{mag}}.gbk"
    threads:
        1
    envmodules:
        {EMAPPER_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    shell:
        """
        emapper2gbk genes -fn {input.fna} -fp {input.faa} -a {input.ann} -o {output.file}
        """

rule m2m:
    input:
        expand(f"{OUTPUT_DIR}/annotating/gbk/{{mag}}/{{mag}}.gbk",mag=mags)
    output:
        sbml=f"{OUTPUT_DIR}/annotating/m2m/sbml/{{mag}}.sbml"
    params:
        indir=f"{OUTPUT_DIR}/annotating/eggnog/",
        outdir=f"{OUTPUT_DIR}/annotating/m2m/"
    threads:
        24
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 100) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 50) * 2 ** (attempt - 1))
    shell:
        """
        module load pathway-tools/27.0 blast/2.13.0 metage2metabo/1.5.3
        m2m recon -g {params.indir} -o {params.outdir} -c {threads}
        """

rule final_network:
    input:
        f"{OUTPUT_DIR}/annotating/m2m/sbml/{{mag}}.sbml"
    output:
        f"{OUTPUT_DIR}/annotating/sbml/{{mag}}.sbml"
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
