####
# Variables parsed from the config.yaml file
####

PACKAGE_DIR = config["package_dir"]

# Software modules
NONPAREIL_MODULE = config["NONPAREIL_MODULE"]

####
# Workflow rules
####

# Estimate metagenomic coverage and diversity from the final metagenomic R1 reads.

rule nonpareil:
    input:
        r1=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/final/{{sample}}_2.fq.gz"
    output:
        npo=f"{OUTPUT_DIR}/preprocessing/nonpareil/{{sample}}.npo",
        stats=f"{OUTPUT_DIR}/preprocessing/nonpareil/{{sample}}_np.tsv"
    params:
        nonpareil_module={NONPAREIL_MODULE},
        prefix=f"{OUTPUT_DIR}/preprocessing/nonpareil/{{sample}}",
        stats_script=f"{PACKAGE_DIR}/workflow/scripts/nonpareil_stats.R"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Estimating Nonpareil coverage and diversity for {wildcards.sample}..."
    shell:
        """
        mkdir -p $(dirname {output.stats:q})
        tmp_dir="${{TMPDIR:-/tmp}}"
        tmp_fastq=$(mktemp "$tmp_dir/{wildcards.sample}.nonpareil.XXXXXX.fq")
        trap 'rm -f "$tmp_fastq"' EXIT

        gunzip -c {input.r1:q} > "$tmp_fastq"
        read_mb=$(python3 -c 'import os, sys; print(os.path.getsize(sys.argv[1]) // 1024 // 1024)' "$tmp_fastq")

        if [ "$read_mb" -gt 150 ]; then
            module load {params.nonpareil_module}
            nonpareil -s "$tmp_fastq" -f fastq -T kmer -t {threads} -b {params.prefix:q}
            Rscript {params.stats_script:q} --sample {wildcards.sample:q} --npo {output.npo:q} --output {output.stats:q}
        else
            : > {output.npo:q}
            printf 'sample\tkappa\tC\tLR\tmodelR\tLRstar\tdiversity\n{wildcards.sample}\t0\t0\t0\t0\t0\t0\n' > {output.stats:q}
        fi
        """
