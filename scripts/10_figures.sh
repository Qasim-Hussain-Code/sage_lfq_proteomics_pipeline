#!/usr/bin/env bash
# Draw every figure in the README from files in results/.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 10_figures.sh [options]

Options:
  --prefix <name>   which arm to plot: sage or maxquant (default: maxquant)
  --force           redraw even if the figures exist
  --help            show this message

Figures whose input file is absent are skipped with a message rather than
failing, so a partial run still produces the plots it can.
USAGE
}

PREFIX="maxquant"; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix) PREFIX="$2"; shift 2 ;;
        --force)  FORCE="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "10_figures: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"
mkdir -p "${FIGURES_DIR}"

if [[ -f "${FIGURES_DIR}/${PREFIX}_volcano_mixed.png" && "${FORCE}" != "true" ]]; then
    echo "10_figures: figures for '${PREFIX}' already drawn, skipping. Use --force."
    exit 0
fi

measure_run "10_figures:${PREFIX}" \
    conda run --no-capture-output --name "${CONDA_ENV_R}" \
    Rscript "${REPO_ROOT}/scripts/figures.R" \
        --results "${RESULTS_DIR}" --figures "${FIGURES_DIR}" --prefix "${PREFIX}" \
    2>&1 | tee "${LOGS_DIR}/10_figures_${PREFIX}.log"

echo "10_figures: figures in ${FIGURES_DIR}"
ls -1 "${FIGURES_DIR}"
