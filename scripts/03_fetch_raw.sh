#!/usr/bin/env bash
# List the deposited files from the PRIDE API, derive the sample sheet from
# the filenames that are actually there, and fetch the RAW files.
#
# The sample sheet is generated, never typed. Filenames in public repositories
# do not always match what a paper says was deposited, and for PXD001584 they
# do not: the authors' own experimental design covers 48 runs across two
# arginine concentrations and two acquisition batches, while only 18 RAW files
# were deposited, all from the 20 micromolar arginine arm. Hardcoding a guess
# at the sample names would have produced a plausible looking sheet describing
# an experiment that is not in the archive.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 03_fetch_raw.sh [options]

Queries the PRIDE API, writes config/samples.tsv, and downloads RAW files.

Options:
  --list-only     write the sample sheet and the manifest, download nothing
  --limit <N>     fetch only the first N runs of the selection (validation runs)
  --subset        one technical replicate per culture (6 runs, not 18)
  --force         re-download files that are already present and correctly sized
  --help          show this message

Checksums: PRIDE returns an empty checksum field for every file in this 2015
submission, so there is no upstream digest to verify against. This script
verifies the byte count against the API's fileSizeBytes and records its own
sha256 in results/raw_manifest.tsv so later runs can detect a changed file.
USAGE
}

LIST_ONLY="false"; LIMIT=""; SUBSET_FLAG=""; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --list-only) LIST_ONLY="true"; shift ;;
        --limit)     LIMIT="$2"; shift 2 ;;
        --subset)    SUBSET_FLAG="true"; shift ;;
        --force)     FORCE="true"; shift ;;
        --help|-h)   usage; exit 0 ;;
        *) echo "03_fetch_raw: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
[[ -n "${SUBSET_FLAG}" ]] && SUBSET="${SUBSET_FLAG}"

RAW_DIR="${DATA_DIR}/raw"
mkdir -p "${RAW_DIR}" "${RESULTS_DIR}" "${CONFIG_DIR}" "${LOGS_DIR}"
SAMPLES="${CONFIG_DIR}/samples.tsv"
MANIFEST="${RESULTS_DIR}/raw_manifest.tsv"
LISTING="${RESULTS_DIR}/pride_file_listing.json"

TMPF=""
# A trap handler that ends on a failed test sets the script exit status
# to 1 even after an explicit "exit 0". Hence the unconditional return.
cleanup() { [[ -n "${TMPF}" && -f "${TMPF}" ]] && rm -f "${TMPF}"; return 0; }
trap cleanup EXIT

# PRIDE API v3. The v2 endpoints were retired; if this 404s, check
# https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html
# for the current path rather than guessing at a URL.
API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects/${PRIDE_ACCESSION}/files?pageSize=200"
if [[ ! -s "${LISTING}" || "${FORCE}" == "true" ]]; then
    echo "03_fetch_raw: querying PRIDE for ${PRIDE_ACCESSION}"
    curl --silent --show-error --fail --location --max-time 300 --output "${LISTING}" "${API}"
else
    echo "03_fetch_raw: reusing cached file listing, delete ${LISTING##*/} to refresh."
fi

# Derive the sample sheet. Filename shape, confirmed against the authors'
# experimentalDesignTemplate.txt in the PRIDE MaxQuant output:
#   20141119_AlCH_<group>_<arginine>_<time>_<culture>_<techrep>_<runid>.raw
# where group 1WT is wild type and 3D8 is the argP transporter deletion,
# arginine is micromolar, culture n3/n4/n5 is the biological replicate and
# techrep 1/2/3 is the injection.
python3 - "${LISTING}" "${SAMPLES}" <<'PY'
import json, re, sys
listing, out = sys.argv[1], sys.argv[2]
files = json.load(open(listing))
pat = re.compile(r'^(\d{8})_AlCH_(\d)(WT|D8)_(\d+)_(\d+h)_n(\d+)_(\d)_(\d+)\.raw$', re.I)
rows = []
for f in files:
    if f.get('fileCategory', {}).get('value') != 'RAW':
        continue
    name = f['fileName']
    m = pat.match(name)
    if not m:
        sys.stderr.write(f"03_fetch_raw: filename does not match expected shape, skipping: {name}\n")
        continue
    date, grpnum, grp, arg, tme, cult, tech, runid = m.groups()
    ftp = next((l['value'] for l in f.get('publicFileLocations', [])
                if 'FTP' in l.get('name', '')), '')
    rows.append(dict(
        run=name[:-4], raw_file=name, sample=f"{grpnum}{grp}_{arg}_{tme}_n{cult}_{tech}",
        genotype=('WT' if grp.upper() == 'WT' else 'argP_KO'),
        arginine_uM=arg, timepoint=tme, biorep=f"n{cult}", techrep=tech,
        run_id=runid, acq_date=date, size_bytes=f.get('fileSizeBytes', 0), url=ftp))
rows.sort(key=lambda r: (r['genotype'] != 'WT', r['biorep'], r['techrep']))
cols = ['run','raw_file','sample','genotype','arginine_uM','timepoint','biorep',
        'techrep','run_id','acq_date','size_bytes','url']
with open(out, 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
tot = sum(int(r['size_bytes']) for r in rows)
print(f"03_fetch_raw: sample sheet has {len(rows)} runs, {tot/2**30:.1f} GiB of RAW")
for g in sorted({r['genotype'] for r in rows}):
    n = sum(1 for r in rows if r['genotype'] == g)
    b = len({r['biorep'] for r in rows if r['genotype'] == g})
    print(f"03_fetch_raw:   {g}: {n} runs across {b} cultures")
PY

echo "03_fetch_raw: wrote ${SAMPLES}"
[[ "${LIST_ONLY}" == "true" ]] && { echo "03_fetch_raw: --list-only, nothing downloaded."; exit 0; }

# Refuse before downloading anything if the selection cannot fit. 00_configure
# projected this from a mean file size; here the exact bytes are known, so the
# check is arithmetic rather than an estimate.
NEEDED_GIB="$(awk -F'\t' -v sub="${SUBSET}" 'NR>1 && (sub!="true" || $8==1) {s+=$11} END {printf "%.1f", s/1073741824}' "${SAMPLES}")"
FREE_GIB="$(df -BG --output=avail "${DATA_DIR}" | tail -1 | tr -dc '0-9')"
echo "03_fetch_raw: selection is ${NEEDED_GIB} GiB of RAW, ${FREE_GIB} GiB free"
echo "03_fetch_raw: RAW is deleted as it converts (04), so peak is one file plus the mzML set"

# Subset selects which runs to fetch; it does not rewrite the sample sheet.
# Keeping technical replicate 1 of each culture is arbitrary. The first
# injection is as defensible as any and it is deterministic.
if [[ "${SUBSET}" == "true" ]]; then
    mapfile -t LINES < <(awk -F'\t' 'NR>1 && $8==1' "${SAMPLES}")
else
    mapfile -t LINES < <(tail -n +2 "${SAMPLES}")
fi
[[ -n "${LIMIT}" ]] && LINES=("${LINES[@]:0:${LIMIT}}")
echo "03_fetch_raw: fetching ${#LINES[@]} file(s)"

: > "${MANIFEST}.tmp"
printf 'run\traw_file\texpected_bytes\tobserved_bytes\tsha256\tstatus\n' > "${MANIFEST}.tmp"

for line in "${LINES[@]}"; do
    IFS=$'\t' read -r run raw_file _ _ _ _ _ _ _ _ size_bytes url <<< "${line}"
    dest="${RAW_DIR}/${raw_file}"
    if [[ -f "${dest}" && "${FORCE}" != "true" ]]; then
        have="$(stat -c%s "${dest}")"
        if [[ "${have}" == "${size_bytes}" ]]; then
            echo "03_fetch_raw: ${raw_file} present and correctly sized, skipping."
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${run}" "${raw_file}" "${size_bytes}" \
                "${have}" "$(sha256sum "${dest}" | cut -d' ' -f1)" "cached" >> "${MANIFEST}.tmp"
            continue
        fi
        echo "03_fetch_raw: ${raw_file} is the wrong size (${have} vs ${size_bytes}), refetching."
    fi
    https_url="${url/ftp:\/\/ftp.pride.ebi.ac.uk/https://ftp.pride.ebi.ac.uk}"
    echo "03_fetch_raw: downloading ${raw_file}"
    TMPF="${dest}.part"
    # aria2c resumes and splits; curl is the fallback so the pipeline does not
    # hard-depend on aria2 being present.
    if command -v aria2c >/dev/null 2>&1; then
        aria2c --quiet --continue=true --max-connection-per-server=4 --split=4 \
               --dir="${RAW_DIR}" --out="${raw_file}.part" "${https_url}"
    else
        curl --silent --show-error --fail --location --continue-at - \
             --output "${TMPF}" "${https_url}"
    fi
    observed="$(stat -c%s "${TMPF}")"
    if [[ "${observed}" != "${size_bytes}" ]]; then
        echo "03_fetch_raw: size mismatch for ${raw_file}: got ${observed}, expected ${size_bytes}" >&2
        rm -f "${TMPF}"; exit 1
    fi
    mv "${TMPF}" "${dest}"; TMPF=""
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${run}" "${raw_file}" "${size_bytes}" \
        "${observed}" "$(sha256sum "${dest}" | cut -d' ' -f1)" "downloaded" >> "${MANIFEST}.tmp"
done

mv "${MANIFEST}.tmp" "${MANIFEST}"
echo "03_fetch_raw: wrote ${MANIFEST}"
