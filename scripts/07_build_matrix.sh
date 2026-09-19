#!/usr/bin/env bash
# Turn Sage's lfq.tsv into the tidy peptide matrix that differential.R reads.
#
# Sage writes one row per quantified precursor and one column per mzML file,
# named after the file. differential.R wants one column per run, named with
# the run names in config/samples.tsv, so that the same script can be pointed
# at the MaxQuant table without modification. The mapping is done here rather
# than in R because it is a string problem, not a statistics problem.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 07_build_matrix.sh [options]

Converts results/sage/lfq.tsv to a tidy peptide matrix plus an annotation table.

Options:
  --lfq <path>      Sage lfq.tsv (default: results/sage/lfq.tsv)
  --output <path>   matrix destination (default: results/peptide_matrix_sage.tsv)
  --force           rebuild even if the matrix exists
  --help            show this message
USAGE
}

LFQ=""; OUTPUT=""; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --lfq)    LFQ="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --force)  FORCE="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "07_build_matrix: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"

LFQ="${LFQ:-${RESULTS_DIR}/sage/lfq.tsv}"
OUTPUT="${OUTPUT:-${RESULTS_DIR}/peptide_matrix_sage.tsv}"
ANNOT="${RESULTS_DIR}/sage_protein_annotation.tsv"
SAMPLES="${CONFIG_DIR}/samples.tsv"

[[ -f "${LFQ}" ]] || die "no ${LFQ}. Run 05_run_sage.sh first (LFQ must be enabled)."
[[ -f "${SAMPLES}" ]] || die "no ${SAMPLES}. Run 03_fetch_raw.sh --list-only."

if [[ -s "${OUTPUT}" && "${FORCE}" != "true" ]]; then
    echo "07_build_matrix: ${OUTPUT##*/} exists, skipping. Use --force to rebuild."
    exit 0
fi

python3 "${REPO_ROOT}/scripts/build_peptide_matrix.py" \
    --source sage --input "${LFQ}" --samples "${SAMPLES}" \
    --output "${OUTPUT}" --annotation "${ANNOT}" \
    --fasta "${DATA_DIR}/fasta/search_database.fasta"

echo "07_build_matrix: wrote ${OUTPUT}"
