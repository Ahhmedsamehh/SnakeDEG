# =====================================================================
#  common.smk  —  shared setup used by every other rule file
# ---------------------------------------------------------------------
#  A .smk file is just Python with extra Snakemake keywords. Anything
#  you can do in Python (import, define functions, build lists) you can
#  do here, and it runs ONCE when Snakemake parses the workflow.
# =====================================================================

import pandas as pd
from snakemake.utils import min_version

# Fail early on very old Snakemake versions instead of with a cryptic error.
min_version("7.0")

# ---------------------------------------------------------------------
#  Load the sample sheet into a pandas DataFrame.
#  `config` is the dict Snakemake built from config/config.yaml.
#  We index by the "sample" column so we can look a row up by name,
#  e.g. samples.loc["control", "fq1"].
# ---------------------------------------------------------------------
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

# A plain list of sample names — used all over the place with expand().
SAMPLES = samples["sample"].tolist()


# ---------------------------------------------------------------------
#  Input functions.
#  A rule's input can be a FUNCTION instead of a string. Snakemake calls
#  it with a `wildcards` object once it knows the concrete sample name,
#  and we return the matching FASTQ path from the sample sheet. This is
#  how each sample gets wired to its own raw reads.
# ---------------------------------------------------------------------
def get_r1(wildcards):
    """Path to the R1 (forward) FASTQ for the requested sample."""
    return samples.loc[wildcards.sample, "fq1"]


def get_r2(wildcards):
    """Path to the R2 (reverse) FASTQ for the requested sample."""
    return samples.loc[wildcards.sample, "fq2"]


# ---------------------------------------------------------------------
#  Constrain the {sample} wildcard to the names we actually know about.
#  Without this, a wildcard could greedily match part of a path and
#  produce confusing "missing rule" errors.
# ---------------------------------------------------------------------
wildcard_constraints:
    sample="|".join(SAMPLES),
