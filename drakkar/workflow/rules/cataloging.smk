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
# Run rules
####

rule assembly:
    input:
        lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_1.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]],
        lambda wildcards: [f"{OUTPUT_DIR}/preprocessing/final/{sample}_2.fq.gz" for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly]]
    output:
        f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    params:
        megahit_module={MEGAHIT_MODULE},
        outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 24) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 200) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    shell:
        """
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
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
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
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
        checkm2_module={CHECKM2_MODULE},
        binette_module={BINETTE_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(sum(preprocessed_mb.get(sample, 1) for sample in assembly_dict_1[wildcards.assembly]) / 1024 * 150) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.checkm2_module} {params.binette_module}
        binette --contig2bin_tables {input.maxbin2} {input.metabat2} --contigs {input.assembly} --outdir {params.outdir}
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
    shell:
        """
        mkdir {params.final}
        mv {params.binette}/final_bins/* {params.final}
        mv {input} {output}
        """
