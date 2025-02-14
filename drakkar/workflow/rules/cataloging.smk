####
# Define config variables
####

MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
METABAT2_MODULE = config["METABAT2_MODULE"]
MAXBIN2_MODULE = config["MAXBIN2_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
BINETTE_MODULE = config["BINETTE_MODULE"]

####
# Calculate file sizes
####

preprocess_mb = calculate_file_sizes(PREPROCESS_DIR)
preprocess_mb = {key.replace('_1.fq.gz', ''): value for key, value in preprocess_mb.items()}
preprocess_mb_total = sum(preprocess_mb.values())

####
# Run rules
####

if "individual" in CATALOGING_MODE:
    rule individual_assembly:
        input:
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 24) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 200) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.megahit_module}
            rm -rf {params.outputdir}
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.outputdir}
            mv {params.outputdir}/final.contigs.fa {output}
            """

    rule individual_assembly_index:
        input:
            f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna"
        output:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.rev.2.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}"
        threads: 1
        resources:
            mem_mb=32*1024,
            runtime=60
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            """

    rule individual_assembly_map:
        input:
            index=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.rev.2.bt2",
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.bam"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 150) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
            """

    rule individual_assembly_map_depth:
        input:
            f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.bam"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.tsv"
        params:
            metabat2_module={METABAT2_MODULE}
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.metabat2_module}
            jgi_summarize_bam_contig_depths --outputDepth {output} {input}
            """

    rule individual_assembly_metabat:
        input:
            assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna",
            depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.tsv"
        output:
            f"{OUTPUT_DIR}/cataloging/metabat2/{{sample}}/bin.1.fa"
        params:
            metabat2_module={METABAT2_MODULE},
            outdir=f"{OUTPUT_DIR}/cataloging/metabat2/{{sample}}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.metabat2_module}
            metabat2 -i {input.assembly} -a {input.depth} -o {params.outdir} -m 1500 --saveCls --noBinOut
            touch {output}
            """

    rule individual_assembly_maxbin:
        input:
            assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna",
            depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{sample}}/{{sample}}.tsv"
        output:
            f"{OUTPUT_DIR}/cataloging/maxbin2/{{sample}}/{{sample}}.summary"
        params:
            maxbin2_module={MAXBIN2_MODULE},
            outdir=f"{OUTPUT_DIR}/cataloging/maxbin2/{{sample}}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.maxbin2_module}
            run_MaxBin.pl -contig {input.assembly} -max_iteration 10 -out {params.outdir} -min_contig_length 1500
            """

    rule individual_binette:
        input:
            metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/{{sample}}/{{sample}}.tsv",
            maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/{{sample}}/{{sample}}.tsv",
            fasta=f"{OUTPUT_DIR}/cataloging/megahit/{{sample}}/{{sample}}.fna"
        output:
            f"{OUTPUT_DIR}/cataloging/binette/{{sample}}/final_bins_quality_reports.tsv"
        params:
            checkm2_module={CHECKM2_MODULE},
            binette_module={BINETTE_MODULE},
            outdir=f"{OUTPUT_DIR}/cataloging/binette/{{sample}}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.checkm2_module} {params.binette_module}
            binette --bin_dirs {input.metabat2} {input.maxbin2} --contigs {input.fasta} --outdir {params.outdir}
            """

if "all" in CATALOGING_MODE:
    rule all_assembly:
        input:
            r1=expand(f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz", sample=samples),
            r2=expand(f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz", sample=samples)
        output:
            f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna"
        params:
            megahit_module={MEGAHIT_MODULE},
            outputdir=f"{OUTPUT_DIR}/cataloging/megahit/all"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(1,int(preprocess_mb_total * 1024 / 40 * 2 ** (attempt - 1))),
            runtime=lambda wildcards, attempt: max(15,int((preprocess_mb_total / 10 * 2 ** (attempt - 1))))
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

    rule all_assembly_index:
        input:
            f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna"
        output:
            index=f"{OUTPUT_DIR}/cataloging/megahit/all/all.rev.2.bt2"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/all/all"
        threads: 1
        resources:
            mem_mb=32*1024,
            runtime=60
        shell:
            """
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
            """

    rule all_assembly_map:
        input:
            index=f"{OUTPUT_DIR}/cataloging/megahit/all/all.rev.2.bt2",
            r1=f"{PREPROCESS_DIR}/{{sample}}_1.fq.gz",
            r2=f"{PREPROCESS_DIR}/{{sample}}_2.fq.gz"
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/all/{{sample}}.bam"
        params:
            bowtie2_module={BOWTIE2_MODULE},
            samtools_module={SAMTOOLS_MODULE},
            basename=f"{OUTPUT_DIR}/cataloging/megahit/all/all"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 150) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.bowtie2_module} {params.samtools_module}
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
            """

    rule all_assembly_map_depth:
        input:
            expand(f"{OUTPUT_DIR}/cataloging/bowtie2/all/{{sample}}.bam", sample=samples)
        output:
            f"{OUTPUT_DIR}/cataloging/bowtie2/all/all.tsv"
        params:
            metabat2_module={METABAT2_MODULE}
        threads: 1
        resources:
            mem_mb=16*1024,
            runtime=60
        shell:
            """
            module load {params.metabat2_module}
            jgi_summarize_bam_contig_depths --outputDepth {output} {input}
            """

    rule all_assembly_metabat:
        input:
            assembly=f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna",
            depth=f"{OUTPUT_DIR}/cataloging/bowtie2/all/all.tsv"
        output:
            f"{OUTPUT_DIR}/cataloging/metabat2/all/bin.1.fa"
        params:
            metabat2_module={METABAT2_MODULE},
            outdir=f"{OUTPUT_DIR}/cataloging/metabat2/all"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.metabat2_module}
            metabat2 -i {input.assembly} -a {input.depth} -o {params.outdir} -m 1500 --saveCls --noBinOut
            touch {output}
            """

    rule all_assembly_maxbin:
        input:
            assembly=f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna",
            depth=f"{OUTPUT_DIR}/cataloging/bowtie2/all/all.tsv"
        output:
            f"{OUTPUT_DIR}/cataloging/maxbin2/all/all.summary"
        params:
            maxbin2_module={MAXBIN2_MODULE},
            outdir=f"{OUTPUT_DIR}/cataloging/maxbin2/all"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocess_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
            runtime=lambda wildcards, attempt: max(10, int(preprocess_mb.get(wildcards.sample, 1) / 1024 * 30) * 2 ** (attempt - 1))
        shell:
            """
            module load {params.maxbin2_module}
            run_MaxBin.pl -contig {input.assembly} -max_iteration 10 -out {params.outdir} -min_contig_length 1500
            """
