# =====================================================================
#  Stage 2 — quality control (FastQC per sample + MultiQC summary)
# ---------------------------------------------------------------------
#  fastp already trims AND reports, but FastQC on the trimmed reads is
#  the community-standard check, and MultiQC then folds every tool's
#  output (fastp, FastQC, HISAT2, featureCounts) into ONE HTML report —
#  the single thing you actually open to judge a run.
# =====================================================================

rule fastqc:
    # We QC the TRIMMED reads. `rules.fastp.output.r1` references the
    # output of the fastp rule by name — cleaner than retyping the path,
    # and it makes the dependency explicit.
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        # FastQC derives output names from the input filename: it strips
        # the .fastq.gz and appends _fastqc.{html,zip}. Because we
        # control the trimmed-read names, these are predictable.
        r1_html="results/qc/fastqc/{sample}_R1.trimmed_fastqc.html",
        r1_zip="results/qc/fastqc/{sample}_R1.trimmed_fastqc.zip",
        r2_html="results/qc/fastqc/{sample}_R2.trimmed_fastqc.html",
        r2_zip="results/qc/fastqc/{sample}_R2.trimmed_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: 2
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir results/qc/fastqc \
            --threads {threads} \
            > {log} 2>&1
        """


rule multiqc:
    # MultiQC must run LAST, so we make it depend on the final QC-bearing
    # outputs of every upstream stage. expand() turns the per-sample
    # template into the full list, e.g. for samples [alzheimers, control]:
    #   results/qc/fastqc/alzheimers_R1.trimmed_fastqc.zip, ...
    input:
        fastqc=expand(
            "results/qc/fastqc/{sample}_R{read}.trimmed_fastqc.zip",
            sample=SAMPLES, read=[1, 2],
        ),
        fastp=expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        hisat2=expand("logs/hisat2/{sample}.log", sample=SAMPLES),
        featurecounts="results/counts/all_samples_counts.tsv.summary",
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        # MultiQC simply scans the folders it is given for recognised log
        # formats. We point it at results/ and logs/ and name the report.
        """
        multiqc results/ logs/ \
            --outdir results/qc \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """
