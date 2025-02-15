rule concatenate_or_link_reads:
    input:
        r1=lambda wildcards: SAMPLE_TO_READS1[wildcards.sample],
        r2=lambda wildcards: SAMPLE_TO_READS2[wildcards.sample]
    output:
        r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
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
