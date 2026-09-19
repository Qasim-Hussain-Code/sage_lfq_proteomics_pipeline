#!/usr/bin/env bash
# Does Sage actually control the false discovery rate it reports?
#
# Build a database of the real proteome plus a paired entrapment proteome,
# search it, and count how many entrapment sequences come back inside the
# reported 1 percent. Any hit to an entrapment sequence is false by
# construction, because those sequences are not in the tube.
#
# The entrapment expansion has to be invisible to the search engine. Sage sees
# 4,200 sequences and has no way to tell the 2,100 shuffled ones from the
# 2,100 real ones: same length, same amino acid composition, same precursor
# masses, same cleavage sites. If it could tell, the whole exercise would be
# circular.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 06_entrapment.sh [options]

Builds the entrapment database, runs the entrapment search, and estimates the
true FDP at PSM, peptide and protein level.

Options:
  --mode <shuffled|foreign>  entrapment construction (default: from project.conf)
  --ratio <r>                entrapment to target ratio (default: 1)
  --thresholds <a,b>         nominal FDR levels to assess (default: 0.01,0.05)
  --build-only               build the database, do not search
  --force                    rebuild and research
  --help                     show this message

The paired estimator of Wen et al. 2025 needs r = 1 and a unique entrapment
partner per target peptide, which only the shuffled mode provides. In foreign
mode the script reports the combined estimator and says so.
USAGE
}

MODE=""; RATIO=""; THRESHOLDS="0.01,0.05"; BUILD_ONLY="false"; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)       MODE="$2"; shift 2 ;;
        --ratio)      RATIO="$2"; shift 2 ;;
        --thresholds) THRESHOLDS="$2"; shift 2 ;;
        --build-only) BUILD_ONLY="true"; shift ;;
        --force)      FORCE="true"; shift ;;
        --help|-h)    usage; exit 0 ;;
        *) echo "06_entrapment: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"

MODE="${MODE:-${ENTRAPMENT_MODE:-shuffled}}"
RATIO="${RATIO:-${ENTRAPMENT_RATIO:-1}}"
# shellcheck disable=SC2153  # DATA_DIR comes from the generated
# project.conf sourced above, which shellcheck cannot follow.
FASTA_DIR="${DATA_DIR}/fasta"
TARGET_FASTA="${FASTA_DIR}/search_database.fasta"
ENT_ONLY="${FASTA_DIR}/entrapment_only.fasta"
COMBINED="${FASTA_DIR}/target_plus_entrapment.fasta"
PAIRING="${FASTA_DIR}/entrapment_pairing.tsv"
OUTDIR="${RESULTS_DIR}/sage_entrapment"

[[ -f "${TARGET_FASTA}" ]] || die "no ${TARGET_FASTA}. Run 02_fetch_fasta.sh."
mkdir -p "${OUTDIR}" "${RESULTS_DIR}" "${LOGS_DIR}"

if [[ -f "${COMBINED}" && "${FORCE}" != "true" ]]; then
    echo "06_entrapment: combined database already built, skipping build."
else
    if [[ "${MODE}" == "shuffled" ]]; then
        echo "06_entrapment: building shuffled entrapment database, seed ${SEED}"
        python3 "${REPO_ROOT}/scripts/build_entrapment_db.py" \
            --target-fasta "${TARGET_FASTA}" \
            --out-fasta "${ENT_ONLY}" \
            --out-pairing "${PAIRING}" \
            --out-stats "${RESULTS_DIR}/entrapment_build_stats.json" \
            --seed "${SEED}"
    else
        FOREIGN="$(find "${FASTA_DIR}" -name '*_foreign.fasta' | head -1)"
        [[ -n "${FOREIGN}" ]] || die "foreign mode needs 02_fetch_fasta.sh --foreign first"
        echo "06_entrapment: using foreign entrapment proteome ${FOREIGN##*/}"
        # No pairing is possible, so equation (4) is unavailable and the
        # assessment falls back to the combined estimator.
        sed 's/^>/>ENT_/' "${FOREIGN}" > "${ENT_ONLY}"
        printf 'target_peptide\tentrapment_peptide\n' > "${PAIRING}"
    fi
    # awk 1 rather than cat: the contaminant FASTA does not end in a
    # newline, and cat would fuse its last sequence line onto the first
    # entrapment header, silently destroying one record and orphaning the
    # sequence after it. This cost me a corrupted database once already.
    awk 1 "${TARGET_FASTA}" "${ENT_ONLY}" > "${COMBINED}"
fi

N_TGT="$(grep -c '^>' "${TARGET_FASTA}")"
N_ENT="$(grep -c '^>' "${ENT_ONLY}")"
N_TOT="$(grep -c '^>' "${COMBINED}")"
echo "06_entrapment: ${N_TGT} target + ${N_ENT} entrapment = ${N_TOT} sequences (mode ${MODE}, r=${RATIO})"

# Record the combined database exactly, so the entrapment result can be tied
# to the sequences that produced it and not merely to the seed.
{
    echo "# Entrapment database provenance"
    echo "mode                    ${MODE}"
    echo "ratio_r                 ${RATIO}"
    echo "seed                    ${SEED}"
    echo "target_sequences        ${N_TGT}"
    echo "entrapment_sequences    ${N_ENT}"
    echo "combined_sequences      ${N_TOT}"
    echo "combined_md5            $(md5sum "${COMBINED}" | cut -d' ' -f1)"
    echo "combined_sha256         $(sha256sum "${COMBINED}" | cut -d' ' -f1)"
    echo "pairing_rows            $(($(wc -l < "${PAIRING}") - 1))"
} > "${RESULTS_DIR}/entrapment_db_provenance.txt"

[[ "${BUILD_ONLY}" == "true" ]] && { echo "06_entrapment: --build-only, stopping here."; exit 0; }

if [[ -s "${OUTDIR}/results.sage.tsv" && "${FORCE}" != "true" ]]; then
    echo "06_entrapment: entrapment search output exists, skipping search."
else
    mapfile -t MZML < <(find "${DATA_DIR}/mzml" -name '*.mzML.gz' -type f | sort)
    [[ "${#MZML[@]}" -gt 0 ]] || die "no mzML found. Run 04_convert_mzml.sh."
    # The entrapment search does not need LFQ: it is asking about
    # identification error, not abundance, and turning quantification off
    # roughly halves the run time.
    ENT_CONFIG="${RESULTS_DIR}/sage_config_entrapment.json"
    python3 - "${CONFIG_DIR}/sage_lfq.json" "${ENT_CONFIG}" "${COMBINED}" "${OUTDIR}" "${MZML[@]}" <<'PY'
import json, sys
base, out, fasta, outdir = sys.argv[1:5]
cfg = json.load(open(base))
cfg['database']['fasta'] = fasta
cfg['output_directory'] = outdir
cfg['mzml_paths'] = sys.argv[5:]
cfg.pop('quant', None)
# predict_rt is only forced on by lfq; with quant removed it can stay off,
# which keeps this search comparable to the main one on identification terms.
cfg['predict_rt'] = False
json.dump(cfg, open(out, 'w'), indent=2)
PY
    echo "06_entrapment: searching ${N_TOT} sequences"
    measure_run "06_entrapment_search" \
        conda run --no-capture-output --name "${CONDA_ENV_MS}" \
        sage --batch-size "${SAGE_BATCH_SIZE:-1}" \
        --disable-telemetry-i-dont-want-to-improve-sage "${ENT_CONFIG}" \
        2>&1 | tee "${LOGS_DIR}/06_entrapment_search.log"
    [[ -f "${OUTDIR}/results.json" ]] && cp "${OUTDIR}/results.json" "${RESULTS_DIR}/sage_entrapment_results.json"
fi

echo "06_entrapment: estimating FDP"
python3 "${REPO_ROOT}/scripts/assess_entrapment.py" \
    --results "${OUTDIR}/results.sage.tsv" \
    --pairing "${PAIRING}" \
    --ratio "${RATIO}" \
    --thresholds "${THRESHOLDS}" \
    --out-json "${RESULTS_DIR}/entrapment_fdr_assessment.json" \
    --out-tsv "${RESULTS_DIR}/entrapment_fdr_assessment.tsv" \
    2>&1 | tee "${LOGS_DIR}/06_entrapment_assessment.log"

echo "06_entrapment: wrote ${RESULTS_DIR}/entrapment_fdr_assessment.tsv"
if grep -q "NOT_CONTROLLED" "${RESULTS_DIR}/entrapment_fdr_assessment.tsv"; then
    echo
    echo "06_entrapment: at least one level is NOT controlled at its nominal rate."
    echo "06_entrapment: this is the headline result, see results/entrapment_fdr_assessment.tsv"
fi
