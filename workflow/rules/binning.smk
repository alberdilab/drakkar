rule binning:
    input:
        "{assembly_dir}/assembly.fasta"
    output:
        "{output_dir}/bins/"
    shell:
        "binning_tool {input} -o {output}"
