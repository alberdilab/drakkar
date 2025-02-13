rule quantification:
    input:
        "{assembly_dir}/assembly.fasta"
    output:
        "{output_dir}/quantification.tsv"
    shell:
        "quant_tool {input} > {output}"
