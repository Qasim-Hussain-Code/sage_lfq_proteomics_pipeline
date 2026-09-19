#!/usr/bin/env bash
# msqrob2 differential abundance on the Sage quantification.
#
# A thin wrapper. All of the statistics live in scripts/differential.R, which
# stage 09 also calls on the MaxQuant table. Keeping one R script rather than
# two is what makes the two results comparable.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 08_differential.sh [options]

Runs scripts/differential.R on the Sage peptide matrix.

Options:
  --matrix <path>       peptide matrix (default: results/peptide_matrix_sage.tsv)
  --prefix <name>       output prefix (default: sage)
  --min-cultures <N>    peptide must be seen in at least N cultures (default: 2)
  --sensitivity         also sweep the culture threshold and report the effect
  --force               rerun even if results exist
  --help                show this message
USAGE
}

MATRIX=""; PREFIX="sage"; MIN_CULTURES="2"; SENS=""; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --matrix)       MATRIX="$2"; shift 2 ;;
        --prefix)       PREFIX="$2"; shift 2 ;;
        --min-cultures) MIN_CULTURES="$2"; shift 2 ;;
        --sensitivity)  SENS="--sensitivity"; shift ;;
        --force)        FORCE="true"; shift ;;
        --help|-h)      usage; exit 0 ;;
        *) echo "08_differential: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"

MATRIX="${MATRIX:-${RESULTS_DIR}/peptide_matrix_sage.tsv}"
[[ -f "${MATRIX}" ]] || die "no ${MATRIX}. Run 07_build_matrix.sh."

OUT="${RESULTS_DIR}/${PREFIX}_protein_results_mixed.tsv"
if [[ -s "${OUT}" && "${FORCE}" != "true" ]]; then
    echo "08_differential: ${OUT##*/} exists, skipping. Use --force to rerun."
    exit 0
fi

measure_run "08_differential:${PREFIX}" \
    conda run --no-capture-output --name "${CONDA_ENV_R}" \
    Rscript "${REPO_ROOT}/scripts/differential.R" \
        --matrix "${MATRIX}" \
        --samples "${CONFIG_DIR}/samples.tsv" \
        --outdir "${RESULTS_DIR}" \
        --prefix "${PREFIX}" \
        --min-cultures "${MIN_CULTURES}" ${SENS:+${SENS}} \
    2>&1 | tee "${LOGS_DIR}/08_differential_${PREFIX}.log"

ANNOT="${RESULTS_DIR}/${PREFIX}_protein_annotation.tsv"
if [[ -f "${ANNOT}" ]]; then
    conda run --no-capture-output --name "${CONDA_ENV_R}" \
        Rscript "${REPO_ROOT}/scripts/biological_readout.R" \
            --results "${OUT}" --annotation "${ANNOT}" \
            --out "${RESULTS_DIR}/${PREFIX}_ribosomal_check.tsv" \
        2>&1 | tee "${LOGS_DIR}/08_ribosomal_${PREFIX}.log"
fi

echo "08_differential: done, results in ${RESULTS_DIR}"
