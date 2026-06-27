# =====================================================================
#  Stage 1 — read trimming & quality filtering with fastp
# ---------------------------------------------------------------------
#  Equivalent to the `fastp ...` block in rnaseq_pipeline.sh, but here it
#  runs once PER SAMPLE automatically, driven by the {sample} wildcard.
# =====================================================================

rule fastp:
    # ---- INPUT ----------------------------------------------------
    # get_r1 / get_r2 (from common.smk) look the raw FASTQ paths up in
    # the sample sheet for whichever {sample} Snakemake is building.
    input:
        r1=get_r1,
        r2=get_r2,
    # ---- OUTPUT ---------------------------------------------------
    # The {sample} in the output paths is the WILDCARD. When Snakemake
    # wants e.g. results/trimmed/control_R1.trimmed.fastq.gz, it sets
    # sample="control" and that value flows into input/log/shell.
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        # fastp's own QC report — later picked up by MultiQC.
        html="results/qc/fastp/{sample}.fastp.html",
        json="results/qc/fastp/{sample}.fastp.json",
    # ---- LOG ------------------------------------------------------
    # Anything written to this path is kept even on success, and is NOT
    # auto-deleted if the job fails (unlike output), so you can debug.
    log:
        "logs/fastp/{sample}.log",
    # ---- THREADS --------------------------------------------------
    # Snakemake passes this number as {threads} and uses it for
    # scheduling (it won't oversubscribe your --cores budget).
    threads: config["threads"]
    # ---- CONDA ----------------------------------------------------
    # Path is relative to THIS .smk file. With --use-conda Snakemake
    # builds the env once and reuses it.
    conda:
        "../envs/fastp.yaml"
    # ---- SHELL ----------------------------------------------------
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            > {log} 2>&1
        """
