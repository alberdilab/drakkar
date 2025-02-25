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
        f"{OUTPUT_DIR}/data/references/{{reference}}.fna"
    message: "Preparing reference genome of {wildcards.reference}..."
    shell:
        """
        if [[ {input} == *.gz ]]; then
            gunzip -c {input} > {output}
        else
            cp {input} {output}
        fi
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
