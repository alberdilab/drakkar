####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
DREP_MODULE = config["DREP_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
COVERM_MODULE = config["COVERM_MODULE"]
MASH_MODULE = config["MASH_MODULE"]
SINGLEM_MODULE = config["SINGLEM_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
DREP_ANI = float(config.get("DREP_ANI", 0.98))
CHECKM2_DB = config["CHECKM2_DB"]

# Annotation databases
SINGLEM_DB = config["SINGLEM_DB"]

####
# Workflow rules
####

rule checkm2_report:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values())
    output:
        report=f"{OUTPUT_DIR}/profiling_genomes/checkm2/quality_report.tsv"
    params:
        checkm2_module={CHECKM2_MODULE},
        checkm2_db={CHECKM2_DB},
        outdir=f"{OUTPUT_DIR}/profiling_genomes/checkm2",
        genome_dir=f"{OUTPUT_DIR}/data/genomes"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Estimating MAG completeness/contamination with CheckM2..."
    shell:
        """
        module load {params.checkm2_module}
        rm -rf {params.outdir}
        rm -rf {params.genome_dir}
        mkdir -p {params.genome_dir}
        for f in {input.genomes}; do
            base=$(basename "$f")
            stem="${{base%.gz}}"
            stem="${{stem%.fa}}"
            stem="${{stem%.fna}}"
            stem="${{stem%.fasta}}"
            out="{params.genome_dir}/${{stem}}.fa"
            if [[ "$f" == *.gz ]]; then
                gunzip -c "$f" > "$out"
            else
                ln -sf "$f" "$out"
            fi
        done
        checkm2 predict --input {params.genome_dir} --output-directory {params.outdir} --threads {threads} --database_path {params.checkm2_db} --extension fa --force
        rm -rf {params.genome_dir}
        """

rule checkm2_metadata:
    input:
        report=f"{OUTPUT_DIR}/profiling_genomes/checkm2/quality_report.tsv"
    output:
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    message: "Formatting CheckM2 report for dRep..."
    run:
        import pandas as pd
        import re
        from pathlib import Path

        report_path = Path(input.report)
        df = pd.read_csv(report_path, sep="\t")
        lower_cols = {c.lower(): c for c in df.columns}
        name_col = lower_cols.get("name") or lower_cols.get("genome") or lower_cols.get("bin")
        comp_col = lower_cols.get("completeness")
        cont_col = lower_cols.get("contamination")
        if not name_col or not comp_col or not cont_col:
            raise ValueError(f"CheckM2 output missing required columns: {list(df.columns)}")

        genome_map = {}
        for genome_path in BINS_TO_FILES.values():
            base = Path(genome_path).name
            base_no_gz = re.sub(r"\.gz$", "", base, flags=re.IGNORECASE)
            base_no_ext = re.sub(r"\.(fa|fna|fasta)$", "", base_no_gz, flags=re.IGNORECASE)
            genome_map.setdefault(base_no_ext, base_no_gz)
            genome_map.setdefault(base_no_gz, base_no_gz)

        genomes = []
        missing = []
        for name in df[name_col].astype(str):
            key = name.strip()
            mapped = genome_map.get(key)
            if not mapped:
                missing.append(key)
                genomes.append(key)
            else:
                genomes.append(mapped)

        if missing:
            raise ValueError(f"CheckM2 names not found in data/genomes: {missing}")

        out_df = pd.DataFrame(
            {
                "genome": genomes,
                "completeness": df[comp_col],
                "contamination": df[cont_col],
            }
        )
        Path(output.metadata).parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(output.metadata, index=False)

checkpoint dereplicate:
    input:
        genomes=expand("{bin_path}", bin_path=BINS_TO_FILES.values()),
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    output:
        Wdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Wdb.csv",
        Cdb=f"{OUTPUT_DIR}/profiling_genomes/drep/data_tables/Cdb.csv"
    params:
        drep_module={DREP_MODULE},
        mash_module={MASH_MODULE},
        outdir=f"{OUTPUT_DIR}/profiling_genomes/drep/",
        ani={DREP_ANI},
        uncompressed_dir=f"{OUTPUT_DIR}/profiling_genomes/drep/uncompressed_genomes"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 10) * 2 ** (attempt - 1))
    message: "Dereplicating bins using dRep..."
    shell:
        """
        module load {params.mash_module} {params.drep_module}
        rm -rf {params.outdir}
        mkdir -p {params.uncompressed_dir}
        files=()
        for f in {input.genomes}; do
            if [[ "$f" == *.gz ]]; then
                out="{params.uncompressed_dir}/$(basename "$f" .gz)"
                gunzip -c "$f" > "$out"
                files+=("$out")
            else
                files+=("$f")
            fi
        done
        dRep dereplicate {params.outdir} -p {threads} -g "${{files[@]}}" -sa {params.ani} --genomeInfo {input.metadata}
        """

# Normalize headers in dereplicated genomes with .fa/.fna/.fasta
rule rename_derep_headers:
    input:
        files=lambda wildcards: get_mag_fna(wildcards)
    output:
        touch(f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes/.headers_renamed")
    localrule: True
    message: "Renaming headers in dereplicated genomes..."
    shell:
        """
        shopt -s nullglob
        for f in {input.files}; do
            base=$(basename "$f")
            genome="${{base%.fa}}"
            genome="${{genome%.fna}}"
            genome="${{genome%.fasta}}"
            awk -v g="$genome" 'BEGIN{{i=0}} /^>/ {{print ">" g "^" i++; next}} {{print}}' "$f" > tmp && mv tmp "$f"
        done
        touch {output}
        """

# Functions to define the input files dynamically.
def get_mag_ids_from_drep(csv_path):
    df = pd.read_csv(csv_path)
    return df["genome"].unique()

def get_mag_fna(wildcards):
    checkpoint_output = checkpoints.dereplicate.get(**wildcards).output[0]
    selected_bins = get_mag_ids_from_drep(checkpoint_output)
    return expand(f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes/{{bin_id}}", bin_id=selected_bins)

rule merge_catalogue:
    input:
        get_mag_fna,
        headers=f"{OUTPUT_DIR}/profiling_genomes/drep/dereplicated_genomes/.headers_renamed"
    output:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    localrule: True
    message: "Merging genomes into a single catalogue..."
    shell:
        """
        cat {input} > {output}
        """

rule index_catalogue:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.fna"
    output:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 1024 * 50) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        """

rule map_to_catalogue:
    input:
        index=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue.rev.1.bt2",
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample]
    output:
        f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=f"{OUTPUT_DIR}/profiling_genomes/catalogue/genome_catalogue"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against genome catalogue..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

rule quantify_reads_catalogue:
    input:
        expand(f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam", sample=samples)
    output:
        f"{OUTPUT_DIR}/profiling_genomes/coverm/coverm.tsv"
    params:
        coverm_module={COVERM_MODULE}
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb / 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 / 5) * 2 ** (attempt - 1))
    message:
        "Generating mapping statistics with..."
    shell:
        """
        module load {params.coverm_module}
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_bases \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """

rule profiling_stats:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.bam"
    output:
        mappedreads=f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.mappedreads",
        mappedbases=f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.mappedbases"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Calculating mapping stats of {wildcards.sample}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        samtools view -b -F12 -@ {threads} {input} | samtools view -c - > {output.mappedreads}
        samtools view -F12 -@ {threads} {input} | awk '{{sum += length($10)}} END {{print sum}}' > {output.mappedbases}
        """

rule profiling_stats_merge:
    input:
        mappedreads=expand(f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.mappedreads", sample=samples),
        mappedbases=expand(f"{OUTPUT_DIR}/profiling_genomes/bowtie2/{{sample}}.mappedbases", sample=samples)
    output:
        f"{OUTPUT_DIR}/profiling_genomes.tsv"
    localrule: True
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=1*1024,
        runtime=5
    message: "Creating profiling genomes stats file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/profiling_genomes_stats.py -r {input.mappedreads} -b {input.mappedbases} -o {output}
        """

rule split_coverm:
    input:
        f"{OUTPUT_DIR}/profiling_genomes/coverm/coverm.tsv"
    output:
        counts=f"{OUTPUT_DIR}/profiling_genomes/final/counts.tsv",
        breadth=f"{OUTPUT_DIR}/profiling_genomes/final/bases.tsv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb = 1*1024,
        runtime = 5
    message: "Generating count and breadth files..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/split_coverm.py {input} {output.counts} {output.breadth}
        """

# Only run if --fraction or -f is used

rule singlem_profile:
    input:
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample]
    output:
        f"{OUTPUT_DIR}/profiling_genomes/singlem/{{sample}}.profile"
    params:
        singlem_module={SINGLEM_MODULE},
        singlem_db={SINGLEM_DB}
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 50) * 2 ** (attempt - 1))
    message: "Running singlem for {wildcards.sample}..."
    shell:
        """
        module load {params.singlem_module}
        singlem pipe --forward {input.r1} --reverse {input.r2} --threads {threads} -p {output} --metapackage {params.singlem_db}
        """

rule singlem_microbial_fraction:
    input:
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample],
        profile=f"{OUTPUT_DIR}/profiling_genomes/singlem/{{sample}}.profile"
    output:
        f"{OUTPUT_DIR}/profiling_genomes/singlem/{{sample}}.tsv"
    params:
        singlem_module={SINGLEM_MODULE},
        singlem_db={SINGLEM_DB}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Running singlem for {wildcards.sample}..."
    shell:
        """
        module load {params.singlem_module}
        singlem microbial_fraction --forward {input.r1} --reverse {input.r2} -p {input.profile} --metapackage {params.singlem_db} > {output}
        """

rule singlem_merge:
    input:
        expand(f"{OUTPUT_DIR}/profiling_genomes/singlem/{{sample}}.tsv", sample=samples)
    output:
        f"{OUTPUT_DIR}/profiling_genomes/singlem/microbial_fraction.tsv"
    localrule: True
    threads: 1
    resources:
        mem_mb=1*1024,
        runtime=5
    message: "Merging singlem outputs..."
    shell:
        """
        awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {input} > {output}
        """
