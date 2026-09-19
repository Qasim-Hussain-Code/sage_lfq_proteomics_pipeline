#!/usr/bin/env bash
# Detect the machine, decide how much of PXD001584 will fit on it, and write
# project.conf. Every later stage sources project.conf instead of re-detecting
# or hardcoding resources, so a run is reproducible from one file.
#
# The disk arithmetic here is the whole point of this stage. The eighteen
# deposited Q Exactive runs are 47.3 GiB of RAW before conversion. A naive
# pipeline downloads all of them, converts all of them, and dies somewhere
# around file fourteen with a full disk and nothing to show. This script
# refuses up front instead.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 00_configure.sh [options]

Writes project.conf in the repository root. Later stages source it.

Options:
  --threads <N>     worker threads for conversion and search (default: nproc)
  --ram <GB>        RAM budget in GiB (default: detected available, not total)
  --disk <GB>       disk budget in GiB (default: detected free on the data volume)
  --dataset <name>  francisella (default) or shen
  --subset          keep one technical replicate per culture: 6 runs, not 18
  --yes             do not prompt, accept the projection and write the file
  --force           overwrite an existing project.conf
  --help            show this message

Notes:
  --subset drops the technical replication and therefore drops the mixed model.
  It is a concession to disk, not a design improvement. See the README.
USAGE
}

THREADS=""
RAM_GB=""
DISK_GB=""
DATASET="francisella"
SUBSET="false"
ASSUME_YES="false"
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --ram)     RAM_GB="$2";  shift 2 ;;
        --disk)    DISK_GB="$2"; shift 2 ;;
        --dataset) DATASET="$2"; shift 2 ;;
        --subset)  SUBSET="true"; shift ;;
        --yes)     ASSUME_YES="true"; shift ;;
        --force)   FORCE="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "00_configure: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONF="${REPO_ROOT}/project.conf"

if [[ -f "${CONF}" && "${FORCE}" != "true" ]]; then
    echo "00_configure: project.conf already exists, skipping. Use --force to rewrite."
    exit 0
fi

case "${DATASET}" in
    francisella|shen) ;;
    *) echo "00_configure: --dataset must be 'francisella' or 'shen'" >&2; exit 2 ;;
esac

# Measured from the PRIDE API on 2026-09-19, not estimated. The eighteen
# deposited runs of PXD001584 total 47.3 GiB, mean 2.63 GiB each. Recording the
# measured mean lets the projection below be arithmetic rather than a guess.
FRANCISELLA_RUNS=18
FRANCISELLA_SUBSET_RUNS=6
FRANCISELLA_MEAN_RAW_GIB="2.63"

# PXD003881 was not enumerated from the API by this script. The run count comes
# from the publication (twenty runs). The mean RAW size is a placeholder that
# 03_fetch_raw.sh replaces with the real listing before anything is downloaded,
# which is why --dataset shen additionally refuses to run without --yes.
SHEN_RUNS=20
SHEN_MEAN_RAW_GIB="3.0"

if [[ "${DATASET}" == "francisella" ]]; then
    N_RUNS=$([[ "${SUBSET}" == "true" ]] && echo "${FRANCISELLA_SUBSET_RUNS}" || echo "${FRANCISELLA_RUNS}")
    MEAN_RAW_GIB="${FRANCISELLA_MEAN_RAW_GIB}"
else
    N_RUNS=$([[ "${SUBSET}" == "true" ]] && echo "$((SHEN_RUNS / 2))" || echo "${SHEN_RUNS}")
    MEAN_RAW_GIB="${SHEN_MEAN_RAW_GIB}"
fi

[[ -n "${THREADS}" ]] || THREADS="$(nproc)"

if [[ -z "${RAM_GB}" ]]; then
    # 'available' rather than 'total'. Total RAM is a lie on a machine that is
    # already running something, and Sage will be killed by the OOM reaper on
    # the basis of what is actually free.
    RAM_GB="$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)"
fi

mkdir -p "${REPO_ROOT}/data"
if [[ -z "${DISK_GB}" ]]; then
    DISK_GB="$(df -BG --output=avail "${REPO_ROOT}/data" | tail -1 | tr -dc '0-9')"
fi

# Conversion ratio. Centroided gzipped mzML from a Q Exactive typically lands
# near a fifth of the RAW size, but "typically" is doing a lot of work in that
# sentence and the real ratio depends on how the file was acquired. This is a
# starting assumption only: 04_convert_mzml.sh measures the true ratio after
# the first file and re-projects against it before converting the second.
MZML_RATIO="0.20"

PROJECTED_MZML_GIB="$(awk -v n="${N_RUNS}" -v m="${MEAN_RAW_GIB}" -v r="${MZML_RATIO}" \
    'BEGIN {printf "%.1f", n*m*r}')"
# Peak disk is all of the retained mzML plus the single RAW being converted at
# that moment. Sage needs every mzML present at once for retention time
# alignment and LFQ, so the mzML total cannot be streamed away.
PROJECTED_PEAK_GIB="$(awk -v z="${PROJECTED_MZML_GIB}" -v m="${MEAN_RAW_GIB}" \
    'BEGIN {printf "%.1f", z+m+2}')"

PROJ_LOG="${REPO_ROOT}/logs/00_configure_projection.txt"
mkdir -p "${REPO_ROOT}/logs"
exec > >(tee -a "${PROJ_LOG}") 2>&1
echo "# ---- $(date -u +%Y-%m-%dT%H:%M:%SZ) ----"
echo "Machine and projection"
echo "  threads available           ${THREADS}"
echo "  RAM available               ${RAM_GB} GiB"
echo "  disk free on data volume    ${DISK_GB} GiB"
echo "  dataset                     ${DATASET} (subset=${SUBSET})"
echo "  runs to process             ${N_RUNS}"
echo "  mean RAW size               ${MEAN_RAW_GIB} GiB"
echo "  projected retained mzML.gz  ${PROJECTED_MZML_GIB} GiB"
echo "  projected peak disk         ${PROJECTED_PEAK_GIB} GiB (mzML + one RAW + 2 GiB headroom)"

SHORTFALL="$(awk -v need="${PROJECTED_PEAK_GIB}" -v have="${DISK_GB}" \
    'BEGIN {d=need-have; printf "%.1f", (d>0 ? d : 0)}')"

if awk -v s="${SHORTFALL}" 'BEGIN {exit !(s > 0)}'; then
    echo
    echo "00_configure: REFUSING to configure a run that cannot finish."
    echo "  projected peak disk ${PROJECTED_PEAK_GIB} GiB exceeds free disk ${DISK_GB} GiB"
    echo "  shortfall           ${SHORTFALL} GiB"
    echo
    if [[ "${SUBSET}" != "true" && "${DATASET}" == "francisella" ]]; then
        echo "  Try --subset, which keeps one technical replicate per culture"
        echo "  (three wild type, three mutant) and needs roughly"
        awk -v m="${FRANCISELLA_MEAN_RAW_GIB}" -v r="${MZML_RATIO}" -v n="${FRANCISELLA_SUBSET_RUNS}" \
            'BEGIN {printf "  %.1f GiB peak instead.\n", n*m*r+m+2}'
        echo "  That drops the technical replication and therefore the mixed model."
    else
        echo "  Free disk, point --disk at a larger volume, or move data/ to one"
        echo "  (data/ may be a symlink to another filesystem; nothing else assumes"
        echo "  it lives inside the repository)."
    fi
    exit 1
fi

if [[ "${ASSUME_YES}" != "true" ]]; then
    read -r -p "Write project.conf with these settings? [y/N] " reply
    case "${reply}" in
        y|Y|yes|YES) ;;
        *) echo "00_configure: aborted, nothing written."; exit 1 ;;
    esac
fi

# One seed for the whole project. The only stochastic step in this pipeline is
# the peptide shuffle that builds the entrapment database (06_entrapment.sh);
# recording the seed here means that database can be rebuilt byte for byte.
SEED=20260919

cat > "${CONF}" <<CONF_EOF
# Written by scripts/00_configure.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)
# Machine-specific. Regenerate with 00_configure.sh --force on a new host.
PROJECT_ROOT="${REPO_ROOT}"
DATA_DIR="${REPO_ROOT}/data"
RESULTS_DIR="${REPO_ROOT}/results"
FIGURES_DIR="${REPO_ROOT}/figures"
LOGS_DIR="${REPO_ROOT}/logs"
CONFIG_DIR="${REPO_ROOT}/config"

THREADS=${THREADS}
RAM_GB=${RAM_GB}
DISK_GB=${DISK_GB}

DATASET="${DATASET}"
SUBSET=${SUBSET}
N_RUNS=${N_RUNS}
MEAN_RAW_GIB=${MEAN_RAW_GIB}
MZML_RATIO=${MZML_RATIO}
PROJECTED_PEAK_GIB=${PROJECTED_PEAK_GIB}

PRIDE_ACCESSION="PXD001584"
UNIPROT_PROTEOME="UP000000762"
UNIPROT_TAXID=401614

# Sage defaults chosen for memory rather than speed. parallel=false and
# batch-size 1 mean Sage holds one file of spectra at a time. On the 1,719
# protein Francisella database the fragment index is small, so the cost of
# this choice is wall clock, not much else. Raise both on a larger machine.
SAGE_BATCH_SIZE=1
SAGE_PARALLEL=false

# Entrapment. Shuffled entrapment at a 1:1 ratio is the default because the
# paired estimator of Wen et al. 2025 requires r=1 with each target peptide
# paired to exactly one entrapment peptide. See 06_entrapment.sh.
ENTRAPMENT_MODE="shuffled"
ENTRAPMENT_RATIO=1
SEED=${SEED}

CONDA_ENV_MS="sage-lfq-ms"
CONDA_ENV_R="sage-lfq-r"
CONF_EOF

echo
echo "00_configure: wrote ${CONF}"
