####
# Define config variables
####

MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
METABAT2_MODULE = config["METABAT2_MODULE"]
MAXBIN2_MODULE = config["MAXBIN2_MODULE"]
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly])) * 5 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 100) * 2 ** (attempt - 1))
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
        mem_mb=32*1024,
        runtime=60
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    message: "Calculating mapping states of assembly {wildcards.assembly}..."
    shell:
        """
        module load {params.metabat2_module}
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """

rule assembly_metabat:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        module load {params.metabat2_module}
        metabat2 -i {input.assembly} -a {input.depth} -o {output} -m 1500 --saveCls --noBinOut
        """

rule assembly_maxbin:
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    run: # Rule run using python to ensure the pipeline does not crush if no bins are generated
        import os
        import subprocess
        try:
            subprocess.run(
                f"module load {params.maxbin2_module} && "
                f"run_MaxBin.pl -contig {input.assembly} -abund {input.depth} "
                f"-max_iteration 10 -out {params.basename} -min_contig_length 1500",
                shell=True,
                check=True
            )
            if os.path.exists(f"{params.basename}.summary"):
                os.rename(f"{params.basename}.summary", output[0])
            else:
                raise FileNotFoundError
        except:
            with open(output[0], "w") as f:
                f.write("")

rule assembly_binette:
    input:
        metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv",
        maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv",
        fasta=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    params:
        diamond_module={DIAMOND_MODULE},
        checkm2_module={CHECKM2_MODULE},
        binette_module={BINETTE_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    message: "Refining bins from assembly {wildcards.assembly} using binette..."
    shell:
        """
        module load {params.diamond_module} {params.checkm2_module} {params.binette_module}

        # Remove empty input files from the list
        VALID_TSV_FILES=()
        for TSV in {input.metabat2} {input.maxbin2}; do
            if [ -s "$TSV" ]; then
                VALID_TSV_FILES+=("$TSV")
            fi
        done

        # Ensure at least one valid TSV file exists
        if [ ${#VALID_TSV_FILES[@]} -eq 0 ]; then
            echo "Error: No valid TSV input files for binette." >&2
            exit 1
        fi

        # Run binette only with non-empty TSV files
        binette --contig2bin_tables "${VALID_TSV_FILES[@]}" \
                --contigs {input.fasta} \
                --outdir {params.outdir} \
                --checkm2_db /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd
        """


rule assembly_final:
    input:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv"
    params:
        binette=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}",
        final=f"{OUTPUT_DIR}/cataloging/final/{{assembly}}"
    threads: 1
    resources:
        mem_mb=8*1024,
        runtime=10
    message: "Outputing final bins from assembly {wildcards.assembly}..."
    shell:
        """
        mkdir {params.final}
        mv {params.binette}/final_bins/* {params.final}
        mv {input} {output}
        """
