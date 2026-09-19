#!/usr/bin/env bash
# Convert Thermo RAW to gzipped centroided mzML, one file at a time, deleting
# each RAW the moment its mzML has been written and verified.
#
# This is the stage that decides whether the pipeline fits on the machine.
# Holding eighteen RAW files and eighteen mzML at once needs about 57 GiB.
# Converting in a stream needs the retained mzML set plus one RAW, which is
# about 14 GiB for the same eighteen runs. The second number fits on a laptop
# and the first does not.
#
# After the first file converts, the real RAW size and the real mzML.gz size
# are known, so the projection for the remaining files becomes arithmetic
# rather than the 0.20 guess in project.conf. If that re-projection exceeds
# free disk the script stops here rather than at file fourteen.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 04_convert_mzml.sh [options]

Streams RAW to gzipped mzML. Fetches each RAW first if it is missing.

Options:
  --limit <N>     convert only the first N runs of the selection
  --runs <a,b>    convert exactly these run names, comma separated. Used for
                  validation runs on a machine that cannot hold the full set
  --subset        one technical replicate per culture (6 runs, not 18)
  --keep-raw      do not delete RAW after conversion (needs far more disk)
  --no-fetch      fail rather than download a RAW that is not already present
  --force         reconvert runs whose mzML already exists
  --help          show this message
USAGE
}

LIMIT=""; SUBSET_FLAG=""; KEEP_RAW="false"; NO_FETCH="false"; FORCE="false"; RUNS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --limit)    LIMIT="$2"; shift 2 ;;
        --runs)     RUNS="$2"; shift 2 ;;
        --subset)   SUBSET_FLAG="true"; shift ;;
        --keep-raw) KEEP_RAW="true"; shift ;;
        --no-fetch) NO_FETCH="true"; shift ;;
        --force)    FORCE="true"; shift ;;
        --help|-h)  usage; exit 0 ;;
        *) echo "04_convert_mzml: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
# shellcheck source=scripts/lib/common.sh
source "${REPO_ROOT}/scripts/lib/common.sh"
[[ -n "${SUBSET_FLAG}" ]] && SUBSET="${SUBSET_FLAG}"

RAW_DIR="${DATA_DIR}/raw"; MZML_DIR="${DATA_DIR}/mzml"
mkdir -p "${RAW_DIR}" "${MZML_DIR}" "${LOGS_DIR}" "${RESULTS_DIR}"
SAMPLES="${CONFIG_DIR}/samples.tsv"
[[ -f "${SAMPLES}" ]] || die "no ${SAMPLES}. Run 03_fetch_raw.sh --list-only first."
CONV_LOG="${RESULTS_DIR}/conversion_sizes.tsv"

PARTIAL=""
# A killed conversion must not leave a truncated mzML.gz behind, because the
# next run would see the file, consider the stage done, and hand Sage a
# corrupt input.
cleanup() { [[ -n "${PARTIAL}" && -f "${PARTIAL}" ]] && rm -f "${PARTIAL}"; return 0; }
trap cleanup EXIT INT TERM

if [[ -n "${RUNS}" ]]; then
    mapfile -t LINES < <(awk -F'\t' -v want=",${RUNS}," 'NR>1 && index(want, ","$1",")' "${SAMPLES}")
    [[ "${#LINES[@]}" -gt 0 ]] || die "none of the requested runs are in the sample sheet"
elif [[ "${SUBSET}" == "true" ]]; then
    mapfile -t LINES < <(awk -F'\t' 'NR>1 && $8==1' "${SAMPLES}")
else
    mapfile -t LINES < <(tail -n +2 "${SAMPLES}")
fi
[[ -n "${LIMIT}" ]] && LINES=("${LINES[@]:0:${LIMIT}}")
TOTAL="${#LINES[@]}"
echo "04_convert_mzml: ${TOTAL} run(s) selected"

[[ -s "${CONV_LOG}" ]] || printf 'run\traw_bytes\tmzml_gz_bytes\tratio\telapsed_s\n' > "${CONV_LOG}"

i=0
for line in "${LINES[@]}"; do
    i=$((i + 1))
    IFS=$'\t' read -r run raw_file _ _ _ _ _ _ _ _ size_bytes _ <<< "${line}"
    mzml="${MZML_DIR}/${run}.mzML.gz"
    raw="${RAW_DIR}/${raw_file}"

    if [[ -f "${mzml}" && "${FORCE}" != "true" ]]; then
        echo "04_convert_mzml: [${i}/${TOTAL}] ${run} already converted, skipping."
        continue
    fi

    if [[ ! -f "${raw}" ]]; then
        [[ "${NO_FETCH}" == "true" ]] && die "${raw_file} missing and --no-fetch was given"
        echo "04_convert_mzml: [${i}/${TOTAL}] fetching ${raw_file}"
        bash "${REPO_ROOT}/scripts/03_fetch_raw.sh" --only "${run}"
    fi

    # Check there is room for this one conversion before starting it.
    need_gib="$(awk -v b="${size_bytes}" -v r="${MZML_RATIO}" 'BEGIN {printf "%.1f", (b*r)/1073741824 + 1}')"
    have_gib="$(free_gib "${DATA_DIR}")"
    if awk -v n="${need_gib}" -v h="${have_gib}" 'BEGIN {exit !(n > h)}'; then
        die "not enough disk to convert ${run}: need about ${need_gib} GiB, have ${have_gib} GiB"
    fi

    echo "04_convert_mzml: [${i}/${TOTAL}] converting ${run}"
    PARTIAL="${mzml}"
    t0="$(date +%s)"
    # -f=1 selects indexed mzML, -g gzips it, -p centroids. Centroided data is
    # what Sage wants and it is several times smaller than profile. The RAW
    # files carry centroided spectra already, so -p is close to a no-op here,
    # but it makes the output independent of how the file was acquired.
    measure_run "04_convert:${run}" \
        conda run --no-capture-output --name "${CONDA_ENV_MS}" \
        ThermoRawFileParser --input="${raw}" --output_file="${mzml}" \
        --format=1 --gzip --noPeakPicking=false >/dev/null
    t1="$(date +%s)"

    # Verify the file before deleting the only other copy of these spectra.
    gzip --test "${mzml}" || die "${mzml##*/} is not valid gzip"
    # An mzML that opens is not necessarily an mzML with spectra in it, so
    # check the root element and that at least one spectrum was written.
    if ! zcat "${mzml}" | head -c 4000 | grep -qE '<(indexedmzML|mzML)'; then
        die "${mzml##*/} does not look like mzML"
    fi
    nspec="$(zcat "${mzml}" | grep -c '<spectrum ' || true)"
    [[ "${nspec}" -gt 0 ]] || die "${mzml##*/} contains no spectra"
    PARTIAL=""

    raw_bytes="$(stat -c%s "${raw}")"
    mzml_bytes="$(stat -c%s "${mzml}")"
    ratio="$(awk -v a="${mzml_bytes}" -v b="${raw_bytes}" 'BEGIN {printf "%.4f", a/b}')"
    printf '%s\t%s\t%s\t%s\t%s\n' "${run}" "${raw_bytes}" "${mzml_bytes}" "${ratio}" "$((t1 - t0))" >> "${CONV_LOG}"
    echo "04_convert_mzml: ${run}: ${nspec} spectra, RAW $((raw_bytes/1048576)) MiB -> mzML.gz $((mzml_bytes/1048576)) MiB (ratio ${ratio})"

    if [[ "${KEEP_RAW}" != "true" ]]; then
        rm -f "${raw}"
        echo "04_convert_mzml: deleted ${raw_file}"
    fi

    # Re-project from the measured ratio, not the configured guess.
    if [[ "${i}" -eq 1 && "${TOTAL}" -gt 1 ]]; then
        remaining=$((TOTAL - 1))
        proj="$(awk -v n="${remaining}" -v m="${mzml_bytes}" -v r="${raw_bytes}" \
            'BEGIN {printf "%.1f", (n*m + r)/1073741824}')"
        have="$(free_gib "${DATA_DIR}")"
        echo "04_convert_mzml: measured ratio ${ratio}; remaining ${remaining} runs project to ${proj} GiB, ${have} GiB free"
        if awk -v p="${proj}" -v h="${have}" 'BEGIN {exit !(p > h)}'; then
            echo "04_convert_mzml: STOPPING after one file." >&2
            echo "04_convert_mzml: the measured conversion ratio is worse than the ${MZML_RATIO} assumed" >&2
            echo "04_convert_mzml: in project.conf, and the remaining runs will not fit." >&2
            echo "04_convert_mzml: shortfall $(awk -v p="${proj}" -v h="${have}" 'BEGIN {printf "%.1f", p-h}') GiB. Try --subset." >&2
            exit 1
        fi
    fi
done

echo "04_convert_mzml: done. mzML retained in ${MZML_DIR}"
du -sh "${MZML_DIR}" 2>/dev/null || true
