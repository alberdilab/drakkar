rule annotation:
    input:
        "{assembly_dir}/assembly.fasta"
    output:
        "{output_dir}/annotations.gff"
    shell:
        "annotation_tool {input} > {output}"
