#!/usr/bin/env bash
# Fetch the search database from UniProt by proteome identifier and record
# enough provenance that the search can be repeated years from now.
#
# A search result without its exact database is not interpretable: peptide to
# protein assignment, the protein level FDR and the razor peptide logic all
# depend on which sequences were present. So this stage writes the UniProt
# release, the download timestamp and both checksums alongside the FASTA, and
# refuses to silently reuse a file it cannot account for.
#
# The original study searched an NCBI-derived FASTA ("NCBI_Fnovicida.fasta",
# from the authors' parameters.txt in the PRIDE submission), not UniProt. That
# is a real difference between this reanalysis and the original and it is
# recorded in the README limitations rather than papered over here.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 02_fetch_fasta.sh [options]

Downloads the target proteome, optionally cRAP contaminants and a foreign
entrapment proteome, and records provenance in results/.

Options:
  --proteome <UPID>   UniProt proteome ID (default: from project.conf)
  --no-contaminants   do not append the cRAP contaminant set
  --foreign           also fetch the foreign entrapment proteome
  --force             re-download even if the FASTA is already present
  --help              show this message
USAGE
}

PROTEOME=""; WITH_CRAP="true"; FETCH_FOREIGN="false"; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --proteome) PROTEOME="$2"; shift 2 ;;
        --no-contaminants) WITH_CRAP="false"; shift ;;
        --foreign)  FETCH_FOREIGN="true"; shift ;;
        --force)    FORCE="true"; shift ;;
        --help|-h)  usage; exit 0 ;;
        *) echo "02_fetch_fasta: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
[[ -n "${PROTEOME}" ]] || PROTEOME="${UNIPROT_PROTEOME}"
[[ "${ENTRAPMENT_MODE:-shuffled}" == "foreign" ]] && FETCH_FOREIGN="true"

# shellcheck disable=SC2153  # DATA_DIR comes from the generated
# project.conf sourced above, which shellcheck cannot follow.
FASTA_DIR="${DATA_DIR}/fasta"
mkdir -p "${FASTA_DIR}" "${RESULTS_DIR}"
TARGET_FASTA="${FASTA_DIR}/${PROTEOME}_target.fasta"
PROV="${RESULTS_DIR}/fasta_provenance.txt"

# Clean up a partial download rather than leaving a truncated FASTA that the
# next stage would happily search against.
# See 03_fetch_raw.sh: a trap ending on a failed test rewrites the exit code.
cleanup() { [[ -n "${TMPF:-}" && -f "${TMPF}" ]] && rm -f "${TMPF}"; return 0; }
trap cleanup EXIT

if [[ -f "${TARGET_FASTA}" && "${FORCE}" != "true" ]]; then
    echo "02_fetch_fasta: ${TARGET_FASTA##*/} already present, skipping download."
else
    echo "02_fetch_fasta: fetching proteome ${PROTEOME} from UniProt"
    TMPF="$(mktemp "${FASTA_DIR}/.dl.XXXXXX")"
    # Canonical sequences only (no isoforms). Bacteria have essentially no
    # isoform annotation, so this changes nothing for Francisella, but it
    # keeps the behaviour explicit if the script is pointed at a eukaryote.
    HDRS="$(mktemp)"
    curl --silent --show-error --fail --location --max-time 600 \
        --dump-header "${HDRS}" \
        --output "${TMPF}" \
        "https://rest.uniprot.org/uniprotkb/stream?query=proteome:${PROTEOME}&format=fasta&includeIsoform=false"
    UNIPROT_RELEASE="$(awk -F': ' 'tolower($1) ~ /x-uniprot-release$/ {gsub(/\r/,"");print $2}' "${HDRS}" | tail -1)"
    UNIPROT_RELEASE_DATE="$(awk -F': ' 'tolower($1) ~ /x-uniprot-release-date$/ {gsub(/\r/,"");print $2}' "${HDRS}" | tail -1)"
    rm -f "${HDRS}"
    [[ -s "${TMPF}" ]] || { echo "02_fetch_fasta: empty download, aborting" >&2; exit 1; }
    mv "${TMPF}" "${TARGET_FASTA}"; TMPF=""
    echo "${UNIPROT_RELEASE:-unknown}" > "${FASTA_DIR}/.uniprot_release"
    echo "${UNIPROT_RELEASE_DATE:-unknown}" > "${FASTA_DIR}/.uniprot_release_date"
fi

N_TARGET="$(grep -c '^>' "${TARGET_FASTA}")"

CRAP_FASTA="${FASTA_DIR}/contaminants.fasta"
# The authors ran MaxQuant with "Include contaminants True", so their
# identifications were competed against a contaminant set. Searching without
# one hands Francisella proteins the spectra that actually came from keratin
# and trypsin, which inflates the identification count and quietly corrupts
# the comparison against their result.
#
# First choice is the universal contaminant library of Frankenfield et al.
# (J Proteome Res 2022, 21(9):2104-2113), whose headers already carry a
# "Cont_" prefix. Second choice is the GPM cRAP set. I reached for cRAP first
# and it failed: on 2026-09-19 ftp.thegpm.org presented a certificate that
# does not match the hostname, so curl refused the connection. Rather than
# pass --insecure to a file that defines what counts as a contaminant, the
# script falls through to the GitHub-hosted library and only then gives up.
CONTAM_PRIMARY="https://raw.githubusercontent.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/main/Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta"
CONTAM_FALLBACK="https://ftp.thegpm.org/fasta/cRAP/crap.fasta"
CONTAM_SOURCE="none"

if [[ "${WITH_CRAP}" == "true" ]]; then
    if [[ -f "${CRAP_FASTA}" && "${FORCE}" != "true" ]]; then
        echo "02_fetch_fasta: contaminants already present, skipping."
        CONTAM_SOURCE="$(cat "${FASTA_DIR}/.contaminant_source" 2>/dev/null || echo cached)"
    else
        echo "02_fetch_fasta: fetching contaminant library"
        TMPF="$(mktemp "${FASTA_DIR}/.dl.XXXXXX")"
        if curl --silent --show-error --fail --location --max-time 300 \
                --output "${TMPF}" "${CONTAM_PRIMARY}"; then
            # Headers already begin with Cont_, so no retagging is needed.
            # 07_build_matrix.sh and differential.R strip on that prefix.
            mv "${TMPF}" "${CRAP_FASTA}"; TMPF=""
            CONTAM_SOURCE="Frankenfield2022_universal"
        elif curl --silent --show-error --fail --location --max-time 300 \
                --output "${TMPF}" "${CONTAM_FALLBACK}"; then
            sed 's/^>/>Cont_/' "${TMPF}" > "${CRAP_FASTA}"
            rm -f "${TMPF}"; TMPF=""
            CONTAM_SOURCE="GPM_cRAP"
        else
            rm -f "${TMPF}"; TMPF=""
            echo "02_fetch_fasta: both contaminant sources failed. Continuing without." >&2
            echo "02_fetch_fasta: the search will be missing a contaminant set, which" >&2
            echo "02_fetch_fasta: inflates identifications. Rerun when network allows." >&2
            WITH_CRAP="false"
        fi
        [[ "${CONTAM_SOURCE}" != "none" ]] && echo "${CONTAM_SOURCE}" > "${FASTA_DIR}/.contaminant_source"
    fi
fi

SEARCH_FASTA="${FASTA_DIR}/search_database.fasta"
if [[ "${WITH_CRAP}" == "true" && -f "${CRAP_FASTA}" ]]; then
    # awk 1, not cat. See 06_entrapment.sh: a source file without a trailing
    # newline makes cat weld two records together.
    awk 1 "${TARGET_FASTA}" "${CRAP_FASTA}" > "${SEARCH_FASTA}"
else
    awk 1 "${TARGET_FASTA}" > "${SEARCH_FASTA}"
fi
N_CRAP=0
[[ -f "${CRAP_FASTA}" && "${WITH_CRAP}" == "true" ]] && N_CRAP="$(grep -c '^>' "${CRAP_FASTA}")"
N_SEARCH="$(grep -c '^>' "${SEARCH_FASTA}")"

if [[ "${FETCH_FOREIGN}" == "true" ]]; then
    # Arabidopsis thaliana. A plant is a defensible foreign entrapment choice
    # for a bacterial sample: it cannot be in the tube, and unlike human it is
    # not itself a routine handling contaminant, so it will not collide with
    # the cRAP set. Homology to Francisella is not zero, which is exactly the
    # pitfall Wen et al. raise about foreign entrapment, and is why the
    # shuffled mode is the default in 06_entrapment.sh.
    FOREIGN_UPID="UP000006548"
    FOREIGN_FASTA="${FASTA_DIR}/${FOREIGN_UPID}_foreign.fasta"
    if [[ -f "${FOREIGN_FASTA}" && "${FORCE}" != "true" ]]; then
        echo "02_fetch_fasta: foreign entrapment proteome already present, skipping."
    else
        echo "02_fetch_fasta: fetching foreign entrapment proteome ${FOREIGN_UPID}"
        TMPF="$(mktemp "${FASTA_DIR}/.dl.XXXXXX")"
        curl --silent --show-error --fail --location --max-time 900 \
            --output "${TMPF}" \
            "https://rest.uniprot.org/uniprotkb/stream?query=proteome:${FOREIGN_UPID}&format=fasta&includeIsoform=false"
        mv "${TMPF}" "${FOREIGN_FASTA}"; TMPF=""
    fi
fi

{
    echo "# Search database provenance"
    echo "# Written by scripts/02_fetch_fasta.sh"
    echo "downloaded_utc          $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "uniprot_proteome        ${PROTEOME}"
    echo "uniprot_taxid           ${UNIPROT_TAXID}"
    echo "uniprot_release         $(cat "${FASTA_DIR}/.uniprot_release" 2>/dev/null || echo unknown)"
    echo "uniprot_release_date    $(cat "${FASTA_DIR}/.uniprot_release_date" 2>/dev/null || echo unknown)"
    echo "target_sequences        ${N_TARGET}"
    echo "contaminants_included   ${WITH_CRAP}"
    echo "contaminant_sequences   ${N_CRAP}"
    echo "contaminant_source      ${CONTAM_SOURCE}"
    echo "contaminant_prefix      Cont_"
    echo "search_db_sequences     ${N_SEARCH}"
    echo "target_fasta_md5        $(md5sum "${TARGET_FASTA}" | cut -d' ' -f1)"
    echo "target_fasta_sha256     $(sha256sum "${TARGET_FASTA}" | cut -d' ' -f1)"
    echo "search_fasta_md5        $(md5sum "${SEARCH_FASTA}" | cut -d' ' -f1)"
    echo "search_fasta_sha256     $(sha256sum "${SEARCH_FASTA}" | cut -d' ' -f1)"
    echo "decoys                  generated internally by Sage, not present in this file"
} > "${PROV}"

echo "02_fetch_fasta: ${N_TARGET} target sequences, ${N_CRAP} contaminants, ${N_SEARCH} in search database"
echo "02_fetch_fasta: wrote ${PROV}"
