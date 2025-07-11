####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
GTDBTK_MODULE = config["GTDBTK_MODULE"]
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]
MMSEQS2_MODULE = config["MMSEQS2_MODULE"]
SIGNALP_MODULE = config["SIGNALP_MODULE"]

# Annotation databases
KEGG_DB = config["KEGG_DB"]
CAZY_DB = config["CAZY_DB"]
AMR_DB = config["AMR_DB"]
PFAM_DB = config["PFAM_DB"]
VFDB_DB = config["VFDB_DB"]
GENOMAD_DB = config["GENOMAD_DB"]

####
# Workflow rules
####

rule prodigal:
    input:
        lambda wildcards: MAGS_TO_FILES[wildcards.mag]
    output:
        nt=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.fna",
        aa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa",
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff"
    params:
        prodigal_module={PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    threads: 1
    message: "Predicting genes of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.prodigal_module}
        if [[ "{input}" == *.gz ]]; then
            gzip -dc {input} | prodigal -d {output.nt} -a {output.aa} -o {output.gff} -f gff
        else
            prodigal -i {input} -d {output.nt} -a {output.aa} -o {output.gff} -f gff
        fi
        """

rule kegg:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        txt=f"{OUTPUT_DIR}/annotating/kegg/{{mag}}.txt",
        tsv=f"{OUTPUT_DIR}/annotating/kegg/{{mag}}.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        db={KEGG_DB}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 60) * 2 ** (attempt - 1))
    threads: 1
    message: "Annotating KEGG orthologs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} -E 1e-10 --noali {params.db} {input}
        """

rule select_kegg:
    input:
        expand(f"{OUTPUT_DIR}/annotating/kegg/{{mag}}.tsv",mag=mags)
    output:
        tsv=f"{OUTPUT_DIR}/annotating/kegg/kegg.tsv",
        csv=f"{OUTPUT_DIR}/annotating/kegg/kegg.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 10 * 2 ** (attempt - 1)
    message: "Selecting and merging all KOs of MAG {wildcards.mag}..."
    shell:
        """
        cat {input} > {output.tsv}
        python {params.package_dir}/workflow/scripts/select_ko.py {output.tsv} {output.csv}
        """

rule cazy:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        txt=f"{OUTPUT_DIR}/annotating/cazy/{{mag}}.txt",
        tsv=f"{OUTPUT_DIR}/annotating/cazy/{{mag}}.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        db={CAZY_DB}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    message: "Annotating CAZYs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input}
        """

rule pfam:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        txt=f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.txt",
        tsv=f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        db={PFAM_DB}
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 30) * 2 ** (attempt - 1))
    message: "Annotating PFAMs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input}
        """

rule vfdb:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        f"{OUTPUT_DIR}/annotating/vfdb/{{mag}}.txt"
    params:
        mmseqs2_module={MMSEQS2_MODULE},
        db={VFDB_DB},
        tmp=f"{OUTPUT_DIR}/annotating/vfdb/{{mag}}"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    message: "Annotating virulence factors of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.mmseqs2_module}
        mmseqs easy-search {input} {params.db} {output} {params.tmp}
        """

rule amr:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        txt=f"{OUTPUT_DIR}/annotating/amr/{{mag}}.txt",
        tsv=f"{OUTPUT_DIR}/annotating/amr/{{mag}}.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        db={AMR_DB}
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    message: "Annotating AMRs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input}
        """

rule signalp:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        f"{OUTPUT_DIR}/annotating/signalp/{{mag}}.txt"
    params:
        signalp_module={SIGNALP_MODULE},
        tmp=f"{OUTPUT_DIR}/annotating/signalp/{{mag}}"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 100) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.signalp_module}
        signalp6 --fastafile {input} --output_dir {params.tmp} --write_procs {threads}
        cat {params.tmp}/output.gff3 | cut -f1,3,6 | awk -F' # |[ \t]+' '!/^#/ {{print $1, $6, $7}}' OFS='\t' > {output}
        """

rule merge_gene_annotations:
    input:
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff",
        kegg=f"{OUTPUT_DIR}/annotating/kegg/{{mag}}.tsv",
        cazy=f"{OUTPUT_DIR}/annotating/cazy/{{mag}}.tsv",
        pfam=f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.tsv",
        vf=f"{OUTPUT_DIR}/annotating/vfdb/{{mag}}.txt",
        amr=f"{OUTPUT_DIR}/annotating/amr/{{mag}}.tsv",
        sp=f"{OUTPUT_DIR}/annotating/signalp/{{mag}}.txt"
    output:
        f"{OUTPUT_DIR}/annotating/final/{{mag}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        kegg=f"{KEGG_DB}.json",
        ec=f"{PFAM_DB}_ec.tsv",
        vf=f"{VFDB_DB}.tsv",
        amr=f"{AMR_DB}.tsv"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    message: "Merging gene annotations of MAG {wildcards.mag}..."
    shell:
        """
        export PATH="/home/jpl786/miniforge3/envs/drakkar/bin:$PATH"
        python {params.package_dir}/workflow/scripts/merge_annotations.py \
            -gff {input.gff} \
            -kegg {input.kegg} \
            -keggdb {params.kegg} \
            -pfam {input.pfam} \
            -ec {params.ec} \
            -cazy {input.cazy} \
            -vf {input.vf} \
            -vfdb {params.vf} \
            -amr {input.amr} \
            -amrdb {params.amr} \
            -signalp {input.sp} \
            -o {output}
        """

rule carbohydrate_clusters:
    input:
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff",
        cazy=f"{OUTPUT_DIR}/annotating/cazy/{{mag}}.tsv",
        pfam=f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.tsv"
    output:
        f"{OUTPUT_DIR}/annotating/cgc/{{mag}}.txt"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 100) * 2 ** (attempt - 1))
    shell:
        """
        python {params.package_dir}/workflow/scripts/detect_cgcs.py \
            -cazy {input.cazy} \
            -pfam {input.pfam}
            -output {output}
        """

rule genomad:
    input:
        lambda wildcards: MAGS_TO_FILES[wildcards.mag]
    output:
        f"{OUTPUT_DIR}/annotating/genomad/{{mag}}/{{mag}}_summary/{{mag}}_virus_summary.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/annotating/genomad/{{mag}}",
        db={GENOMAD_DB}
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 100) * 2 ** (attempt - 1))
    shell:
        """
        output/GPB_bin_000118_summary/GPB_bin_000118_virus_summary.tsv
        module load genomad/1.11.0
        genomad end-to-end \
            -t {threads} \
            {input} \
            {params.outdir} \
            {params.db}
        """

rule final_gene_annotation_table:
    input:
        expand(f"{OUTPUT_DIR}/annotating/final/{{mag}}.tsv",mag=mags)
    output:
        f"{OUTPUT_DIR}/annotating/gene_annotations.tsv.xz"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    message: "Generating final gene annotation file..."
    shell:
        """
        cat {input} | xz -c > {output}
        """
