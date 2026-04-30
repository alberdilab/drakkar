####
# Config variables
####

PACKAGE_DIR = config["package_dir"]
project_name = config["project_name"]

####
# Workflow rules
####

rule create_report:
    output:
        f"{OUTPUT_DIR}/drakkar_report.html"
    localrule: True
    params:
        package_dir={PACKAGE_DIR},
        project_name={project_name}
    threads: 1
    resources:
        mem_mb=1*1024,
        runtime=5
    message: "Initialising DRAKKAR report..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/create_report.py -p {params.project_name} -o {output}
        """

rule concatenate_or_link_reads:
    input:
        r1=lambda wildcards: SAMPLE_TO_READS1[wildcards.sample],
        r2=lambda wildcards: SAMPLE_TO_READS2[wildcards.sample]
    output:
        r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
    message: "Preparing input of {wildcards.sample}..."
    shell:
        """
        if [ $(echo {input.r1} | wc -w) -gt 1 ]; then
            cat {input.r1} > {output.r1};
        else
            ln -s $(realpath {input.r1}) {output.r1};
        fi

        if [ $(echo {input.r2} | wc -w) -gt 1 ]; then
            cat {input.r2} > {output.r2};
        else
            ln -s $(realpath {input.r2}) {output.r2};
        fi
        """

rule prepare_reference:
    input:
        lambda wildcards: REFERENCE_TO_FILE[wildcards.reference]
    output:
        fna=f"{OUTPUT_DIR}/data/references/{{reference}}.fna",
        ready=f"{OUTPUT_DIR}/data/references/{{reference}}.index.ready"
    params:
        bowtie2_module=config["BOWTIE2_MODULE"],
        extract_script=f"{PACKAGE_DIR}/workflow/scripts/extract_reference_index.py",
        reference_dir=f"{OUTPUT_DIR}/data/references",
        basename=f"{OUTPUT_DIR}/data/references/{{reference}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Preparing reference genome of {wildcards.reference}..."
    shell:
        """
        reference_input={input:q}

        case "$reference_input" in
            *.tar|*.tar.gz|*.tgz|*.tar.bz2|*.tbz2|*.tar.xz|*.txz)
                python {params.extract_script:q} \
                    --archive "$reference_input" \
                    --reference {wildcards.reference:q} \
                    --output-dir {params.reference_dir:q}
                ;;
            *.gz)
                gunzip -c "$reference_input" > {output.fna:q}
                module load {params.bowtie2_module}
                bowtie2-build {output.fna:q} {params.basename:q}
                ;;
            *)
                cp "$reference_input" {output.fna:q}
                module load {params.bowtie2_module}
                bowtie2-build {output.fna:q} {params.basename:q}
                ;;
        esac
        if [ ! -s {output.fna:q} ]; then
            echo "ERROR: Reference preparation did not create the required FASTA file." >&2
            exit 1
        fi
        if [ ! -s {params.basename:q}.rev.1.bt2 ] && [ ! -s {params.basename:q}.rev.1.bt2l ]; then
            echo "ERROR: Reference preparation did not create required FASTA and Bowtie2 index files." >&2
            exit 1
        fi
        touch {output.ready:q}
        """

rule concatenate_or_link_preprocessed:
    input:
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample]
    output:
        r1=f"{OUTPUT_DIR}/preprocessed/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessed/final/{{sample}}_2.fq.gz"
    message: "Preparing input of {wildcards.sample}..."
    shell:
        """
        if [ $(echo {input.r1} | wc -w) -gt 1 ]; then
            cat {input.r1} > {output.r1};
        else
            ln -s $(realpath {input.r1}) {output.r1};
        fi

        if [ $(echo {input.r2} | wc -w) -gt 1 ]; then
            cat {input.r2} > {output.r2};
        else
            ln -s $(realpath {input.r2}) {output.r2};
        fi
        """
