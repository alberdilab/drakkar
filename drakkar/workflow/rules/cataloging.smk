####
# Define config variables
####

PACKAGE_DIR = config["package_dir"]
MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
METABAT2_MODULE = config["METABAT2_MODULE"]
MAXBIN2_MODULE = config["MAXBIN2_MODULE"]
SEMIBIN2_MODULE = config["SEMIBIN2_MODULE"]
DIAMOND_MODULE = config["DIAMOND_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
BINETTE_MODULE = config["BINETTE_MODULE"]

####
# Run rules
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
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb) * 2 ** (attempt - 1))
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
        r1=lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_1.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]],
        r2=lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_2.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]]
    output:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{assembly}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=lambda wildcards: f"{OUTPUT_DIR}/cataloging/megahit/{wildcards.assembly}/{wildcards.assembly}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping reads to assembly {wildcards.assembly}..."
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule assembly_map_depth:
    input:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{assembly}}.bam"
    output:
        metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.depth",
        maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.depth"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
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
        depth=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.depth"
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
        depth=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv"
    params:
        maxbin2_module={MAXBIN2_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    shell:
        """
        module load {params.maxbin2_module}
        run_MaxBin.pl -contig {input.assembly} -abund {input.depth} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        """

#not active currently
rule semibin2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/{{assembly}}.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/recluster_bins_info.tsv"
    params:
        semibin2_module={SEMIBIN2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        module load {params.semibin2_module}
        SemiBin2 single_easy_bin -i {input.assembly} --depth-metabat2 {input.depth} -o {output} -m 1500
        """

checkpoint assembly_binette:
    input:
        metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv",
        maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv",
        fasta=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    params:
        diamond_module=DIAMOND_MODULE,
        checkm2_module=CHECKM2_MODULE,
        binette_module=BINETTE_MODULE,
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 40) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 4) * 2 ** (attempt - 1))
    message: "Refining bins from assembly {wildcards.assembly} using binette..."
    shell:
        """
        module load {params.diamond_module} {params.checkm2_module} {params.binette_module}

        # Define input files
        METABAT2="{input.metabat2}"
        MAXBIN2="{input.maxbin2}"

        # Remove empty input files from the list
        VALID_TSV_FILES=""
        if [ -s "$METABAT2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $METABAT2"
        fi
        if [ -s "$MAXBIN2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $MAXBIN2"
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
                --checkm2_db /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd
        """

# Functions to define the input files dynamically.
def get_bin_ids_from_tsv(tsv_path):                                                                     #
    df = pd.read_csv(tsv_path, sep="\t")                                                                #
    return df["bin_id"].unique()                                                                        #

def get_bin_fna_sep(wildcards):                                                                         # all this is probably unnecessary
    checkpoint_output = checkpoints.assembly_binette.get(**wildcards).output[0]                         # and the rename_bins input can be
    cluster_ids = get_bin_ids_from_tsv(checkpoint_output)                                               # simplified, as bin_id is created in the main Snakemake
    return f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins/bin_{wildcards.bin_id}.fa"

rule rename_bins:
    input:
        bin=lambda wildcards: get_bin_fna_sep(wildcards)
    output:
        bin=f"{OUTPUT_DIR}/cataloging/final/{{assembly}}/{{assembly}}_bin_{{bin_id}}.fa"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(2, int(input.size_mb) * 2 ** (attempt - 1))
    message: "Renaming bin {wildcards.bin_id} from assembly {wildcards.assembly}..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/rename_bins.py {wildcards.assembly} {input.bin} {output.bin}
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
