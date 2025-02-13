rule assembly:
    input:
        "{input_dir}/reads.fq"
    output:
        "{output_dir}/assembly.fasta"
    shell:
        "assembly_tool {input} > {output}"
