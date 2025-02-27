####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
GTDBTK_MODULE = config["GTDBTK_MODULE"]
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]

# Annotation databases
KEGG_DB = config["KEGG_DB"]
CAZY_DB = config["CAZY_DB"]
AMR_DB = config["AMR_DB"]
PFAM_DB = config["PFAM_DB"]
VFDB_DB = config["VFDB_DB"]

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
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    threads: 1
    message: "Predicting genes of MAG {wildards.mag}..."
    shell:
        """
        module load {params.prodigal_module}
        prodigal -i {input} -d {output.nt} -a {output.aa} -o {output.gff} -f gff
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
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
        f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.txt"
    params:
        mmseqs2_module={MMSEQS2_MODULE},
        db={VFDB_DB},
        tmp=f"{OUTPUT_DIR}/annotating/pfam/{wildcards.mag}"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Annotating virulence factors of MAG {wildcards.mag}..."
    shell:
        """
        module load m{params.mmseqs2_module}
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Annotating AMRs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input}
        """

rule merge_annotations:
    input:
        gff=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.gff",
        kegg=f"{OUTPUT_DIR}/annotating/kegg/{{mag}}.tsv",
        cazy=f"{OUTPUT_DIR}/annotating/cazy/{{mag}}.tsv",
        pfam=f"{OUTPUT_DIR}/annotating/pfam/{{mag}}.tsv",
        vfdb=f"{OUTPUT_DIR}/annotating/vfdb/{{mag}}.txt",
        amrdb=f"{OUTPUT_DIR}/annotating/amr/{{mag}}.tsv"
        #sp="results/signalp/{genome}.txt"
    output:
        f"{OUTPUT_DIR}/annotating/final/{{mag}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        kegg=f"{KEGG_DB}.json",
        ec=f"{KEGG_DB}_ec.tsv",
        vf=f"{VFDB_DB}.tsv",
        amr=f"{AMR_DB}.tsv"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Merging annotations of MAG {wildcards.mag}..."
    shell:
        """
        python {params.package_dir}/scripts/merge_annotations.py \
            -gff {input.gff} \
            -keggdb {params.kegg} \
            -kegg {input.kegg} \
            -pfam {input.pfam} \
            -ec {params.ec} \
            -cazy {input.cazy} \
            -vf {input.vfdb} \
            -vfdb {params.vf} \
            -amr {input.amr} \
            -amrdb {params.amr} \
            -o {output}
        """

rule final_annotation_table:
    input:
        expand(f"{OUTPUT_DIR}/annotating/final/{{mag}}.tsv",mag=mags)
    output:
        f"{OUTPUT_DIR}/annotating/gene_annotations.tsv.xz"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    message: "Generating final gene annotation file..."
    shell:
        """
        cat {input.tsv_files} | xz -c > {output}
        """
