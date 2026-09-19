#!/usr/bin/env bash
# Run the Sage search with label free quantification over every converted run.
#
# Sage needs all of the mzML present at once. Label free quantification aligns
# retention times across runs and traces MS1 peaks between them, so there is
# no streaming version of this stage: whatever the conversion stage retained
# is what gets searched, together.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 05_run_sage.sh [options]

Searches data/mzml/*.mzML.gz against data/fasta/search_database.fasta.

Options:
  --fasta <path>        override the search database
  --output <dir>        override the output directory (default: results/sage)
  --config <path>       override the base config (default: config/sage_lfq.json)
  --batch-size <N>      files loaded and searched at once (default: 1)
  --parallel            set parallel:true in the config, faster and hungrier
  --allow-telemetry     do not pass --disable-telemetry to Sage
  --force               rerun even if results.sage.tsv already exists
  --help                show this message

Telemetry: Sage sends a small startup report unless told not to. This script
passes --disable-telemetry by default. See the README for what is sent.
USAGE
}

FASTA_OVERRIDE=""; OUT_OVERRIDE=""; CONFIG_OVERRIDE=""
BATCH_OVERRIDE=""; PARALLEL_FLAG="false"; ALLOW_TELEMETRY="false"; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)      FASTA_OVERRIDE="$2"; shift 2 ;;
        --output)     OUT_OVERRIDE="$2"; shift 2 ;;
        --config)     CONFIG_OVERRIDE="$2"; shift 2 ;;
        --batch-size) BATCH_OVERRIDE="$2"; shift 2 ;;
        --parallel)   PARALLEL_FLAG="true"; shift ;;
        --allow-telemetry) ALLOW_TELEMETRY="true"; shift ;;
        --force)      FORCE="true"; shift ;;
        --help|-h)    usage; exit 0 ;;
        *) echo "05_run_sage: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"

MZML_DIR="${DATA_DIR}/mzml"
FASTA="${FASTA_OVERRIDE:-${DATA_DIR}/fasta/search_database.fasta}"
OUTDIR="${OUT_OVERRIDE:-${RESULTS_DIR}/sage}"
BASE_CONFIG="${CONFIG_OVERRIDE:-${CONFIG_DIR}/sage_lfq.json}"
BATCH_SIZE="${BATCH_OVERRIDE:-${SAGE_BATCH_SIZE:-1}}"

[[ -f "${FASTA}" ]] || die "no search database at ${FASTA}. Run 02_fetch_fasta.sh."
[[ -f "${BASE_CONFIG}" ]] || die "no config at ${BASE_CONFIG}"
mkdir -p "${OUTDIR}" "${LOGS_DIR}"

mapfile -t MZML < <(find "${MZML_DIR}" -name '*.mzML.gz' -type f | sort)
[[ "${#MZML[@]}" -gt 0 ]] || die "no mzML in ${MZML_DIR}. Run 04_convert_mzml.sh."

if [[ -s "${OUTDIR}/results.sage.tsv" && "${FORCE}" != "true" ]]; then
    echo "05_run_sage: ${OUTDIR}/results.sage.tsv exists, skipping. Use --force to rerun."
    exit 0
fi

# Sage drops any FASTA entry whose accession starts with decoy_tag when it is
# generating decoys internally. A database that already contains "rev_"
# accessions would lose them silently, so check rather than hope.
DECOY_TAG="$(python3 -c "import json;print(json.load(open('${BASE_CONFIG}'))['database'].get('decoy_tag','rev_'))")"
COLLISIONS="$(grep -c "^>${DECOY_TAG}" "${FASTA}" || true)"
if [[ "${COLLISIONS}" -gt 0 ]]; then
    die "${COLLISIONS} FASTA entries already begin with the decoy tag '${DECOY_TAG}'. Sage would discard them. Change decoy_tag."
fi

# Write the config that actually runs, with paths resolved. The base config in
# config/ stays readable and path-independent; this copy is the record.
USED_CONFIG="${RESULTS_DIR}/sage_config_used.json"
python3 - "${BASE_CONFIG}" "${USED_CONFIG}" "${FASTA}" "${OUTDIR}" "${PARALLEL_FLAG}" "${MZML[@]}" <<'PY'
import json, sys
base, out, fasta, outdir, parallel = sys.argv[1:6]
mzml = sys.argv[6:]
cfg = json.load(open(base))
cfg['database']['fasta'] = fasta
cfg['output_directory'] = outdir
cfg['mzml_paths'] = mzml
if parallel == 'true':
    cfg['parallel'] = True
json.dump(cfg, open(out, 'w'), indent=2)
print(f"05_run_sage: {len(mzml)} mzML file(s), database {fasta}")
PY

SAGE_ARGS=(--batch-size "${BATCH_SIZE}")
if [[ "${ALLOW_TELEMETRY}" != "true" ]]; then
    SAGE_ARGS+=(--disable-telemetry)
fi

echo "05_run_sage: searching with batch-size ${BATCH_SIZE}, parallel=${PARALLEL_FLAG}"
measure_run "05_sage_search" \
    conda run --no-capture-output --name "${CONDA_ENV_MS}" \
    sage "${SAGE_ARGS[@]}" "${USED_CONFIG}" \
    2>&1 | tee "${LOGS_DIR}/05_sage_search.log"

# Sage writes results.json recording every parameter it resolved, including
# the defaults you did not set. That file is the machine-readable methods
# section, so it is kept rather than left in the output directory to be
# overwritten by the entrapment search.
if [[ -f "${OUTDIR}/results.json" ]]; then
    cp "${OUTDIR}/results.json" "${RESULTS_DIR}/sage_results.json"
    echo "05_run_sage: kept ${RESULTS_DIR}/sage_results.json"
fi

for f in results.sage.tsv lfq.tsv; do
    [[ -f "${OUTDIR}/${f}" ]] || echo "05_run_sage: warning, expected ${f} was not written" >&2
done

if [[ -f "${OUTDIR}/results.sage.tsv" ]]; then
    n_psm="$(($(wc -l < "${OUTDIR}/results.sage.tsv") - 1))"
    echo "05_run_sage: ${n_psm} PSM rows written"
fi
