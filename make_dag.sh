#!/usr/bin/env bash
# Regenerate the workflow DAG figure (images/dag.png + images/dag.pdf) with a
# clean, professional style.
#
# It asks Snakemake for the rule graph as Graphviz DOT, restyles it (uniform
# palette, Helvetica, filled rounded boxes), and renders it with `dot`.
#
# Requires graphviz:  conda install -c conda-forge graphviz   (or apt-get install graphviz)
# Usage (from the personal_pipeline/ directory):  bash make_dag.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

if ! command -v dot >/dev/null 2>&1; then
    echo "ERROR: graphviz 'dot' not found. Install it with:" >&2
    echo "  conda install -c conda-forge graphviz   # or: apt-get install -y graphviz" >&2
    exit 1
fi

mkdir -p images

# --rulegraph = one node per rule (clean). Use --dag for one node per job.
restyle() {
    snakemake --rulegraph 2>/dev/null \
    | sed -E \
        -e 's/color = "[0-9. ]*"/color="#34495E"/g' \
        -e 's/style="rounded"/style="rounded,filled", fillcolor="#D6EAF8"/g' \
        -e 's/node\[[^]]*\];/node[shape=box, fontname="Helvetica", fontsize=14, margin="0.3,0.16", penwidth=2];/' \
        -e 's/edge\[[^]]*\];/edge[color="#7F8C8D", penwidth=1.6, arrowsize=0.9];/'
}

restyle | dot -Tpng -Gdpi=200 -Grankdir=TB -Granksep=0.6 -Gnodesep=0.5 -o images/dag.png
restyle | dot -Tpdf            -Grankdir=TB -Granksep=0.6 -Gnodesep=0.5 -o images/dag.pdf

echo "Wrote images/dag.png and images/dag.pdf"
