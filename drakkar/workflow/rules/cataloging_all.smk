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

rule all_assembly:
    input:
        r1=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz", sample=samples),
        r2=expand(f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz", sample=samples)
    output:
        f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna"
    params:
        megahit_module={MEGAHIT_MODULE},
        outputdir=f"{OUTPUT_DIR}/cataloging/megahit/all"
    threads: 8
    resources:
        mem_mb=max(8 * 1024,int(preprocessed_mb_total * 1024 / 40 * 2 ** (attempt - 1))),
        runtime=max(60,int((preprocessed_mb_total / 10 * 2 ** (attempt - 1))))
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
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
    output:
        f"{OUTPUT_DIR}/cataloging/bowtie2/all/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/megahit/all/all"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: max(8*1024, int(preprocessed_mb.get(wildcards.sample, 1) * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: max(10, int(preprocessed_mb.get(wildcards.sample, 1) / 1024 * 150) * 2 ** (attempt - 1))
    shell:
        """
        module load {params.bowtie2_module} {params.samtools_module}
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule all_assembly_map_depth:
    input:
        expand(f"{OUTPUT_DIR}/cataloging/bowtie2/all/{{sample}}.bam", sample=samples)
    output:
        metabat2=f"{OUTPUT_DIR}/cataloging/metabat2/all/all.depth",
        maxbin2=f"{OUTPUT_DIR}/cataloging/maxbin2/all/all.depth"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=16*1024,
        runtime=60
    shell:
        """
        module load {params.metabat2_module}
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """

rule all_assembly_metabat:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna",
        depth=f"{OUTPUT_DIR}/cataloging/metabat2/all/all.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/metabat2/all/all.tsv"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=24*1024,
        runtime=120
    shell:
        """
        module load {params.metabat2_module}
        metabat2 -i {input.assembly} -a {input.depth} -o {output} -m 1500 --saveCls --noBinOut
        """

rule all_assembly_maxbin:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/all/all.fna",
        depth=f"{OUTPUT_DIR}/cataloging/maxbin2/all/all.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/all/all.tsv"
    params:
        maxbin2_module={MAXBIN2_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/maxbin2/all/all"
    threads: 1
    resources:
        mem_mb=24*1024,
        runtime=120
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
