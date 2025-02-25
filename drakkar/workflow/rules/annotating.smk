####
# Define config variables
####

PACKAGE_DIR = config["package_dir"]
GTDBTK_MODULE = config["GTDBTK_MODULE"]
PRODIGAL_MODULE = config["PRODIGAL_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]

rule gtdbtk_input:
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    params:
        package_dir={PACKAGE_DIR},
        names=expand("{mag_name}", mag_name=MAGS_TO_FILES.keys()),
        paths=expand("{mag_path}", mag_path=MAGS_TO_FILES.values())
    localrule: True
    threads: 1
    resources:
        mem_mb = 1*1024,
        runtime = 5
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/gtdbtk_input.py --names {params.names} --paths {params.paths} --output {output}
        """

rule gtdbtk:
    input:
        f"{OUTPUT_DIR}/annotating/gtdbtk/mag_input.tsv"
    output:
        f"{OUTPUT_DIR}/annotating/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        gtdbtk_module={GTDBTK_MODULE},
        data=f"/datasets/globe_databases/gtdbtk_db/20241001",
        outdir=f"{OUTPUT_DIR}/profiling_genomes/gtdbtk/",
        tmpdir=f"{OUTPUT_DIR}/profiling_genomes/tmp/"
    threads: 24
    resources:
        mem_mb=lambda wildcards, input, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, input, attempt: 120 * 2 ** (attempt - 1)
    message: "Annotating taxonomy using GTDBTK..."
    shell:
        """
        module load {params.gtdbtk_module}
        export GTDBTK_DATA_PATH={params.data}
        mkdir -p {params.tmpdir}
        gtdbtk classify_wf \
            --batchfile {input} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --skip_ani_screen
        """

rule prodigal:
    input:
        lambda wildcards: MAGS_TO_FILES[wildcards.mag]
    output:
        nt=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.fna",
        aa=f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    params:
        prodigal_module={PRODIGAL_MODULE}
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024) * 2 ** (attempt - 1))
    threads: 1
    message: "Predicting genes using Prodigal..."
    shell:
        """
        module load {params.prodigal_module}
        prodigal -i {input} -d {output.nt} -a {output.aa}
        """

rule kofams:
    input:
        f"{OUTPUT_DIR}/annotating/prodigal/{{mag}}.faa"
    output:
        txt=f"{OUTPUT_DIR}/annotating/kofams/{{mag}}.txt",
        tsv=f"{OUTPUT_DIR}/annotating/kofams/{{mag}}.tsv"
    params:
        hmmer_module={HMMER_MODULE},
        database="/maps/datasets/globe_databases/dram/20240606/kofam_profiles.hmm"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 5) * 2 ** (attempt - 1))
    threads: 1
    message: "Annotating KEGG orthologs of MAG {wildcards.mag}..."
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} -E 1e-10 --noali {params.database} {input}
        """

rule select_kegg:
    input:
        expand(f"{OUTPUT_DIR}/annotating/kofams/{{mag}}.tsv",mag=mags)
    output:
        tsv=f"{OUTPUT_DIR}/annotating/kofams/kegg.tsv",
        csv=f"{OUTPUT_DIR}/annotating/kofams/kegg.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 10 * 2 ** (attempt - 1)
    message: "Selecting and merging all KOs..."
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
        db="/maps/datasets/globe_databases/dram/20240606/dbCAN-HMMdb-V11.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 10 * 2 ** (attempt - 1)
    shell:
        """
        module load {params.hmmer_module}
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input}
        """
