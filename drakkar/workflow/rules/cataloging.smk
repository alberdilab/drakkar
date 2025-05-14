####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
METABAT2_MODULE = config["METABAT2_MODULE"]
MAXBIN2_MODULE = config["MAXBIN2_MODULE"]
BEDTOOLS_MODULE = config["BEDTOOLS_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]
SEMIBIN2_MODULE = config["SEMIBIN2_MODULE"]
DIAMOND_MODULE = config["DIAMOND_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
BINETTE_MODULE = config["BINETTE_MODULE"]

# Databases
CHECKM2_DB = config["CHECKM2_DB"]

####
# Workflow rules
####

rule assembly:
    input:
        r1=lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_1.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]],
        r2=lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_2.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]]
    output:
        f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    params:
        megahit_module={MEGAHIT_MODULE},
        outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1020*1024,max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Assembling {wildcards.assembly}..."
    shell:
        """
        module load {params.megahit_module}
        rm -rf {params.outputdir}

        # Convert input list to a comma-separated string
        R1_FILES=$(echo {input.r1} | tr ' ' ',')
        R2_FILES=$(echo {input.r2} | tr ' ' ',')

        megahit \
            -t {threads} \
            --verbose \
            --min-contig-len 1500 \
            -1 $R1_FILES -2 $R2_FILES \
            -o {params.outputdir}
        mv {params.outputdir}/final.contigs.fa {output}
        """

rule assembly_index:
    input:
        f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        index=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.rev.2.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Indexing assembly {wildcards.assembly}..."
    shell:
        """
        module load {params.bowtie2_module}
        bowtie2-build {input} {params.basename}
        """

rule assembly_map:
    input:
        index=lambda wildcards: f"{OUTPUT_DIR}/cataloging/megahit/{wildcards.assembly}/{wildcards.assembly}.rev.2.bt2",
        r1=lambda wildcards: f"{OUTPUT_DIR}/preprocessing/final/{wildcards.sample}_1.fq.gz",
        r2=lambda wildcards: f"{OUTPUT_DIR}/preprocessing/final/{wildcards.sample}_2.fq.gz"
    output:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=lambda wildcards: f"{OUTPUT_DIR}/cataloging/megahit/{wildcards.assembly}/{wildcards.assembly}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping reads to assembly {wildcards.assembly}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule assembly_map_depth:
    input:
        lambda wildcards: [
            f"{OUTPUT_DIR}/cataloging/bowtie2/{wildcards.assembly}/{sample}.bam"
            for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]
        ]
    output:
        metabat2=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_metabat.depth",
        maxbin2=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_maxbin.depth"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Calculating mapping states of assembly {wildcards.assembly}..."
    shell:
        """
        module load {params.metabat2_module}
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """

rule metabat2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        module load {params.metabat2_module}
        metabat2 -i {input.assembly} -a {input.depth} -o {output} -m 1500 --saveCls --noBinOut
        """

rule maxbin2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_maxbin.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.summary"
    params:
        maxbin2_module={MAXBIN2_MODULE},
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    shell:
        """
        module load {params.maxbin2_module} {params.bowtie2_module}
        rm -rf {params.basename}*
        run_MaxBin.pl -contig {input.assembly} -abund {input.depth} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        """

rule maxbin2_table:
    input:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.summary"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        fastadir=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    shell:
        """
        python {params.package_dir}/workflow/scripts/fastas_to_bintable.py -d {params.fastadir} -e fasta -o {output}
        """

rule semibin2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        bam=lambda wildcards: [
            f"{OUTPUT_DIR}/cataloging/bowtie2/{wildcards.assembly}/{sample}.bam"
            for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]
            ]
    output:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/contig_bins.tsv"
    params:
        semibin2_module={SEMIBIN2_MODULE},
        hmmer_module={HMMER_MODULE},
        bedtools_module={BEDTOOLS_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 50) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using semibin2..."
    shell:
        """
        module load {params.semibin2_module} {params.bedtools_module} {params.hmmer_module}
        SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o {params.outdir} -m 1500 -t {threads} --compression none
        """

rule semibin2_table:
    input:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/contig_bins.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        fastadir=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/output_bins"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(5, int(input.size_mb / 5) * 2 ** (attempt - 1))
    shell:
        """
        tail -n +2 {input} > {output}
        """

checkpoint binette:
    input:
        metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv",
        maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv",
        semibin2=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/{{assembly}}.tsv",
        fasta=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    params:
        checkm_db = {CHECKM2_DB},
        diamond_module = {DIAMOND_MODULE},
        checkm2_module = {CHECKM2_MODULE},
        binette_module = {BINETTE_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 100) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 2) * 2 ** (attempt - 1))
    message: "Refining bins from assembly {wildcards.assembly} using binette..."
    shell:
        """
        module load {params.diamond_module} {params.checkm2_module} {params.binette_module}

        # Define input files
        METABAT2="{input.metabat2}"
        MAXBIN2="{input.maxbin2}"
        SEMIBIN2="{input.semibin2}"

        # Remove empty input files from the list
        VALID_TSV_FILES=""
        if [ -s "$METABAT2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $METABAT2"
        fi
        if [ -s "$MAXBIN2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $MAXBIN2"
        fi
        if [ -s "$SEMIBIN2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $SEMIBIN2"
        fi

        # Ensure at least one valid TSV file exists
        if [ -z "$VALID_TSV_FILES" ]; then
            echo "Error: No valid TSV input files for binette." >&2
            exit 1
        fi

        # Run binette only with non-empty TSV files
        binette --contig2bin_tables $VALID_TSV_FILES \
                --contigs {input.fasta} \
                --outdir {params.outdir} \
                --checkm2_db {params.checkm_db}
        """

# Regenerate the bin_id wildcard based on the checkpoint results
def get_bin_fna_sep(wildcards):
    checkpoint_output = checkpoints.binette.get(**wildcards).output[0]
    cluster_ids = get_bin_ids_from_tsv(checkpoint_output)
    return f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins/bin_{wildcards.bin_id}.fa"

rule rename_bins:
    input:
        lambda wildcards: get_bin_fna_sep(wildcards)
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}/{{assembly}}_bin_{{bin_id}}.fa"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(2, int(input.size_mb) * 2 ** (attempt - 1))
    message: "Renaming bin {wildcards.bin_id} from assembly {wildcards.assembly}..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/rename_bins.py {wildcards.assembly} {input} {output}
        """

rule move_metadata:
    input:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv"
    threads: 1
    resources:
        mem_mb=8*1024,
        runtime=10
    message: "Exporting bin metadata from assembly {wildcards.assembly}..."
    shell:
        """
        cp {input} {output}
        """

rule all_bins:
    input:
        expand(f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv", assembly=assemblies)
    output:
        paths=f"{OUTPUT_DIR}/cataloging/final/all_bin_paths.txt",
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(2, int(input.size_mb) * 2 ** (attempt - 1))
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/all_bin_paths.py {input} -o {output.paths}
        python {params.package_dir}/workflow/scripts/all_bin_metadata.py {input} -o {output.metadata}
        """
