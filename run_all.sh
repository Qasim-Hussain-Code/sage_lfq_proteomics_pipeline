#!/usr/bin/env bash
# Run the pipeline end to end, or from any stage onwards.
#
# Every stage is idempotent: re-running one that has already finished prints
# what it is skipping and returns. So --from is a convenience, not a
# requirement, and an interrupted run is restarted by running this again.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

usage() {
    cat <<'USAGE'
Usage: run_all.sh [options]

Stages, in order:
  00_configure      detect threads, RAM and disk; project the footprint; write project.conf
  01_install        conda environments for the MS tools and for R/Bioconductor
  02_fetch_fasta    UniProt proteome by proteome ID, plus the contaminant library
  03_fetch_raw      query the PRIDE API, write config/samples.tsv, fetch RAW
  04_convert_mzml   RAW to gzipped mzML one at a time, deleting each RAW as it goes
  05_run_sage       Sage search with label free quantification
  06_entrapment     entrapment database, entrapment search, FDR assessment
  07_build_matrix   Sage lfq.tsv to a tidy peptide matrix
  08_differential   msqrob2 modelling of the Sage quantification
  09_compare        the authors' MaxQuant table through the same workflow, then compare
  10_figures        every figure in the README

Options:
  --from <stage>       start at this stage (name or number, e.g. 07 or 07_build_matrix)
  --only <stage>       run exactly one stage
  --dataset <name>     francisella (default) or shen
  --subset             one technical replicate per culture: 6 runs, not 18
  --limit <N>          process only the first N runs (validation runs)
  --threads <N>        override detected thread count
  --disk <GB>          override detected free disk
  --maxquant-only      run only the stages that need no RAW files:
                       00, 01, 09 (MaxQuant arm), 10. Needs about 1 GB of disk.
  --skip-install       do not run 01_install
  --yes                do not prompt in 00_configure
  --force              pass --force to each stage, redoing completed work
  --help               show this message

Examples:
  ./run_all.sh --yes                       full pipeline, eighteen runs
  ./run_all.sh --subset --yes              six runs, no mixed model
  ./run_all.sh --maxquant-only --yes       no RAW files, reproduces the MaxQuant arm
  ./run_all.sh --from 08                   re-model without re-searching
USAGE
}

FROM=""; ONLY=""; DATASET="francisella"; SUBSET=""; LIMIT=""
THREADS=""; DISK=""; MQ_ONLY="false"; SKIP_INSTALL="false"; ASSUME_YES=""; FORCE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)          FROM="$2"; shift 2 ;;
        --only)          ONLY="$2"; shift 2 ;;
        --dataset)       DATASET="$2"; shift 2 ;;
        --subset)        SUBSET="--subset"; shift ;;
        --limit)         LIMIT="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --disk)          DISK="$2"; shift 2 ;;
        --maxquant-only) MQ_ONLY="true"; shift ;;
        --skip-install)  SKIP_INSTALL="true"; shift ;;
        --yes)           ASSUME_YES="--yes"; shift ;;
        --force)         FORCE="--force"; shift ;;
        --help|-h)       usage; exit 0 ;;
        *) echo "run_all: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

STAGES=(00_configure 01_install 02_fetch_fasta 03_fetch_raw 04_convert_mzml
        05_run_sage 06_entrapment 07_build_matrix 08_differential 09_compare 10_figures)

normalise() {
    local want="$1" s
    for s in "${STAGES[@]}"; do
        [[ "${s}" == "${want}" || "${s%%_*}" == "${want}" ]] && { echo "${s}"; return 0; }
    done
    echo "run_all: no such stage '${want}'" >&2
    echo "run_all: stages are ${STAGES[*]}" >&2
    return 1
}

START_IDX=0
if [[ -n "${FROM}" ]]; then
    FROM="$(normalise "${FROM}")"
    for i in "${!STAGES[@]}"; do [[ "${STAGES[$i]}" == "${FROM}" ]] && START_IDX="${i}"; done
fi
if [[ -n "${ONLY}" ]]; then ONLY="$(normalise "${ONLY}")"; fi

if [[ "${MQ_ONLY}" == "true" ]]; then
    RUN_SET=(00_configure 01_install 09_compare 10_figures)
else
    RUN_SET=("${STAGES[@]:${START_IDX}}")
fi
[[ -n "${ONLY}" ]] && RUN_SET=("${ONLY}")

banner() { printf '\n==== %s ====\n' "$1"; }

for stage in "${RUN_SET[@]}"; do
    case "${stage}" in
        00_configure)
            banner "00_configure"
            # shellcheck disable=SC2086
            bash scripts/00_configure.sh --dataset "${DATASET}" ${SUBSET} ${ASSUME_YES} \
                ${THREADS:+--threads ${THREADS}} ${DISK:+--disk ${DISK}} ${FORCE}
            ;;
        01_install)
            [[ "${SKIP_INSTALL}" == "true" ]] && { echo "run_all: skipping 01_install"; continue; }
            banner "01_install"
            bash scripts/01_install.sh
            ;;
        02_fetch_fasta)
            banner "02_fetch_fasta"
            # shellcheck disable=SC2086
            bash scripts/02_fetch_fasta.sh ${FORCE}
            ;;
        03_fetch_raw)
            banner "03_fetch_raw"
            # The sample sheet is always written. Fetching is left to 04,
            # which interleaves it with conversion so only one RAW is ever
            # on disk at a time.
            # shellcheck disable=SC2086
            bash scripts/03_fetch_raw.sh --list-only ${SUBSET}
            ;;
        04_convert_mzml)
            banner "04_convert_mzml"
            # shellcheck disable=SC2086
            bash scripts/04_convert_mzml.sh ${SUBSET} ${LIMIT:+--limit ${LIMIT}} ${FORCE}
            ;;
        05_run_sage)
            banner "05_run_sage"
            # shellcheck disable=SC2086
            bash scripts/05_run_sage.sh ${FORCE}
            ;;
        06_entrapment)
            banner "06_entrapment"
            # shellcheck disable=SC2086
            bash scripts/06_entrapment.sh ${FORCE}
            ;;
        07_build_matrix)
            banner "07_build_matrix"
            # shellcheck disable=SC2086
            bash scripts/07_build_matrix.sh ${FORCE}
            ;;
        08_differential)
            banner "08_differential"
            # shellcheck disable=SC2086
            bash scripts/08_differential.sh --sensitivity ${FORCE}
            ;;
        09_compare)
            banner "09_compare_maxquant"
            # shellcheck disable=SC2086
            if [[ "${MQ_ONLY}" == "true" ]]; then
                bash scripts/09_compare_maxquant.sh --maxquant-only --sensitivity ${FORCE}
            else
                bash scripts/09_compare_maxquant.sh --sensitivity ${FORCE}
            fi
            ;;
        10_figures)
            banner "10_figures"
            # shellcheck disable=SC2086
            bash scripts/10_figures.sh --prefix maxquant ${FORCE}
            if [[ -s results/sage_protein_results_mixed.tsv ]]; then
                # shellcheck disable=SC2086
                bash scripts/10_figures.sh --prefix sage ${FORCE}
            fi
            ;;
    esac
done

banner "done"
if [[ -s logs/stage_metrics.tsv ]]; then
    echo "Measured cost per stage (logs/stage_metrics.tsv):"
    column -t -s$'\t' logs/stage_metrics.tsv
fi
