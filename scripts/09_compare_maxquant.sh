#!/usr/bin/env bash
# Run the identical msqrob2 workflow on the authors' own MaxQuant table and
# compare it to the Sage result.
#
# This is the part of the project that is worth publishing. The two arms share
# their spectra and their statistics and differ only in the search engine and
# the quantification, so a disagreement localises to those. The MaxQuant arm
# also stands on its own: it needs no RAW files and no search, only a 10 MB
# table, which means it runs on any machine in about a minute.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 09_compare_maxquant.sh [options]

Fetches the authors' MaxQuant peptides.txt, runs differential.R on it, and
compares the result against the Sage arm.

Options:
  --maxquant-only   run the MaxQuant arm, skip the comparison (no Sage needed)
  --min-cultures N  peptide must be seen in at least N cultures (default: 2)
  --sensitivity     sweep the culture threshold in the MaxQuant arm
  --force           refetch and rerun
  --help            show this message

The peptide table is the authors' MaxQuant output restricted to the eighteen
deposited runs, redistributed by statOmics with the accessions re-annotated
from NCBI gi identifiers to RefSeq WP_ identifiers. It is fetched, not
vendored, so this repository does not redistribute their data.
USAGE
}

MQ_ONLY="false"; MIN_CULTURES="2"; SENS=""; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --maxquant-only) MQ_ONLY="true"; shift ;;
        --min-cultures)  MIN_CULTURES="$2"; shift 2 ;;
        --sensitivity)   SENS="--sensitivity"; shift ;;
        --force)         FORCE="true"; shift ;;
        --help|-h)       usage; exit 0 ;;
        *) echo "09_compare_maxquant: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"

MQ_DIR="${DATA_DIR}/maxquant"
mkdir -p "${MQ_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}"
MQ_PEPTIDES="${MQ_DIR}/peptides.txt"
MQ_MATRIX="${MQ_DIR}/peptide_matrix_maxquant.tsv"
MQ_URL="https://raw.githubusercontent.com/statOmics/msqrob2data/refs/heads/main/dda/francisellaFull/peptides.txt"

TMPF=""
cleanup() { [[ -n "${TMPF}" && -f "${TMPF}" ]] && rm -f "${TMPF}"; return 0; }
trap cleanup EXIT

if [[ -s "${MQ_PEPTIDES}" && "${FORCE}" != "true" ]]; then
    echo "09_compare_maxquant: peptides.txt already present, skipping download."
else
    echo "09_compare_maxquant: fetching the authors' MaxQuant peptide table"
    TMPF="$(mktemp "${MQ_DIR}/.dl.XXXXXX")"
    curl --silent --show-error --fail --location --max-time 600 --output "${TMPF}" "${MQ_URL}"
    mv "${TMPF}" "${MQ_PEPTIDES}"; TMPF=""
fi
{
    echo "source_url            ${MQ_URL}"
    echo "downloaded_utc        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "sha256                $(sha256sum "${MQ_PEPTIDES}" | cut -d' ' -f1)"
    echo "rows                  $(($(wc -l < "${MQ_PEPTIDES}") - 1))"
} > "${RESULTS_DIR}/maxquant_source_provenance.txt"

if [[ -s "${MQ_MATRIX}" && "${FORCE}" != "true" ]]; then
    echo "09_compare_maxquant: MaxQuant matrix present, skipping rebuild."
else
    python3 "${REPO_ROOT}/scripts/build_peptide_matrix.py" \
        --source maxquant --input "${MQ_PEPTIDES}" --samples "${CONFIG_DIR}/samples.tsv" \
        --output "${MQ_MATRIX}" --annotation "${RESULTS_DIR}/maxquant_protein_annotation.tsv"
fi

if [[ -s "${RESULTS_DIR}/maxquant_protein_results_mixed.tsv" && "${FORCE}" != "true" ]]; then
    echo "09_compare_maxquant: MaxQuant differential results present, skipping."
else
    measure_run "09_differential:maxquant" \
        conda run --no-capture-output --name "${CONDA_ENV_R}" \
        Rscript "${REPO_ROOT}/scripts/differential.R" \
            --matrix "${MQ_MATRIX}" --samples "${CONFIG_DIR}/samples.tsv" \
            --outdir "${RESULTS_DIR}" --prefix maxquant \
            --min-cultures "${MIN_CULTURES}" ${SENS:+${SENS}} \
        2>&1 | tee "${LOGS_DIR}/09_differential_maxquant.log"
    conda run --no-capture-output --name "${CONDA_ENV_R}" \
        Rscript "${REPO_ROOT}/scripts/biological_readout.R" \
            --results "${RESULTS_DIR}/maxquant_protein_results_mixed.tsv" \
            --annotation "${RESULTS_DIR}/maxquant_protein_annotation.tsv" \
            --out "${RESULTS_DIR}/maxquant_ribosomal_check.tsv" \
        2>&1 | tee "${LOGS_DIR}/09_ribosomal_maxquant.log"
fi

[[ "${MQ_ONLY}" == "true" ]] && { echo "09_compare_maxquant: --maxquant-only, stopping."; exit 0; }

SAGE_RES="${RESULTS_DIR}/sage_protein_results_mixed.tsv"
if [[ ! -s "${SAGE_RES}" ]]; then
    echo "09_compare_maxquant: no Sage results at ${SAGE_RES##*/}."
    echo "09_compare_maxquant: the MaxQuant arm is complete and its outputs are in results/."
    echo "09_compare_maxquant: run 05 to 08 to produce the Sage arm, then rerun this stage."
    exit 0
fi

# UniProt accession to RefSeq, so the two engines' protein identifiers can be
# compared. Sage reports UniProt because that is the database it searched;
# the authors reported RefSeq because that is what they searched.
MAPPING="${RESULTS_DIR}/uniprot_to_refseq.tsv"
if [[ -s "${MAPPING}" && "${FORCE}" != "true" ]]; then
    echo "09_compare_maxquant: accession mapping present, skipping fetch."
else
    echo "09_compare_maxquant: fetching UniProt to RefSeq cross-references"
    TMPF="$(mktemp)"
    curl --silent --show-error --fail --location --max-time 300 --output "${TMPF}" \
        "https://rest.uniprot.org/uniprotkb/stream?query=proteome:${UNIPROT_PROTEOME}&format=tsv&fields=accession,xref_refseq"
    # Explode the semicolon-delimited RefSeq list into one row per pair and
    # drop the version suffix, which the two sources do not agree on.
    awk -F'\t' 'NR>1 {
        n = split($2, a, ";")
        for (i = 1; i <= n; i++) {
            gsub(/^ +| +$/, "", a[i]); sub(/\.[0-9]+$/, "", a[i])
            if (a[i] != "") print $1 "\t" a[i]
        }
    }' "${TMPF}" | sort -u | awk 'BEGIN {print "uniprot\trefseq"} {print}' > "${MAPPING}"
    rm -f "${TMPF}"; TMPF=""
fi
echo "09_compare_maxquant: $(($(wc -l < "${MAPPING}") - 1)) UniProt to RefSeq pairs"

conda run --no-capture-output --name "${CONDA_ENV_R}" \
    Rscript "${REPO_ROOT}/scripts/compare_engines.R" \
        --sage "${SAGE_RES}" \
        --maxquant "${RESULTS_DIR}/maxquant_protein_results_mixed.tsv" \
        --mapping "${MAPPING}" \
        --sage-matrix "${RESULTS_DIR}/peptide_matrix_sage.tsv" \
        --maxquant-matrix "${MQ_MATRIX}" \
        --outdir "${RESULTS_DIR}" \
    2>&1 | tee "${LOGS_DIR}/09_compare_engines.log"

echo "09_compare_maxquant: wrote ${RESULTS_DIR}/engine_comparison.tsv"
