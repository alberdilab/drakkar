####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
SEQKIT_MODULE = config["SEQKIT_MODULE"]

# seqkit pairs mates by hashing the record ID verbatim, and takes the leading
# non-space characters as that ID. Raw BGI/DNBSEQ headers carry no space and
# end in "/1" and "/2" (@V300026712L2C001R0010000372/1), so the suffix lands
# inside the ID and no mate ever matches, silently yielding zero pairs. The
# regexp strips it. Illumina headers separate the mate field with a space and
# need no such handling, so the flag is only passed for BGI.
PLATFORM = str(config.get("platform", "illumina")).lower()
PAIR_ID_REGEXP = r"--id-regexp '^(\S+)\/[12]'" if PLATFORM == "bgi" else ""

####
# Workflow rules
####

# Sanitize each read file independently with seqkit sana to remove malformed records,
# then re-pair the two files so only reads present in both are retained.
# Outputs sanitized paired FASTQ files and a stats file recording the retained read count.

rule seqkit_sana:
    input:
        r1=f"{OUTPUT_DIR}/data/reads/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/data/reads/{{sample}}_2.fq.gz"
    output:
        r1=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_1.fq.gz",
        r2=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_2.fq.gz",
        stats=f"{OUTPUT_DIR}/preprocessing/seqkit/{{sample}}_sana.tsv"
    params:
        seqkit_module=SEQKIT_MODULE,
        pair_id_regexp=PAIR_ID_REGEXP,
        workdir=lambda wildcards: f"{OUTPUT_DIR}/preprocessing/seqkit/.tmp_{wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 5)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 512)) * 2 ** (attempt - 1))
    message: "Sanitizing FASTQ files for sample {wildcards.sample}..."
    shell:
        """
        module purge
        module load {params.seqkit_module}
        mkdir -p {params.workdir}/paired
        seqkit sana {input.r1} -o {params.workdir}/sana_1.fq.gz
        seqkit sana {input.r2} -o {params.workdir}/sana_2.fq.gz
        seqkit pair {params.pair_id_regexp} \
            -1 {params.workdir}/sana_1.fq.gz \
            -2 {params.workdir}/sana_2.fq.gz \
            -O {params.workdir}/paired
        mv {params.workdir}/paired/sana_1.fq.gz {output.r1}
        mv {params.workdir}/paired/sana_2.fq.gz {output.r2}
        count=$(zcat {output.r1} | awk 'NR%4==1{{count++}} END{{print count * 2 + 0}}')
        if [ "$count" -eq 0 ]; then
            echo "ERROR: seqkit pair matched no read pairs for {wildcards.sample}. The mates share no read IDs, which usually means the FASTQ headers carry a mate suffix that seqkit folds into the ID. If these are BGI/DNBSEQ reads, rerun with --platform bgi." >&2
            rm -rf {params.workdir}
            exit 1
        fi
        printf "sample\tseqkit_sana_reads\n{wildcards.sample}\t$count" > {output.stats}
        rm -rf {params.workdir}
        """
