#!/usr/bin/env bash
# The pre-publication checks, run as one command.
#
# Point 5 of the brief is that every number in the README has to be traceable
# to a file in results/ or logs/. This script prints the numbers alongside
# the file each one came from, so that claim can be checked rather than
# asserted. It also runs the mechanical checks: shellcheck, banned prose, em
# dashes, emojis, and the size of anything tracked by git.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"
fail=0
note() { printf '\n== %s ==\n' "$1"; }

note "1. shellcheck"
if command -v shellcheck >/dev/null 2>&1; then
    if shellcheck --shell=bash run_all.sh scripts/*.sh scripts/lib/*.sh; then
        echo "clean, no findings at any severity"
    else
        echo "shellcheck reported findings"; fail=1
    fi
else
    echo "shellcheck not installed, skipped"
fi

note "2. every script responds to --help"
for s in run_all.sh scripts/[0-9]*.sh; do
    if bash "$s" --help >/dev/null 2>&1; then printf '  ok   %s\n' "$s"
    else printf '  FAIL %s\n' "$s"; fail=1; fi
done

note "3. em dashes and emojis in tracked text"
# U+2014 em dash, U+2013 en dash, and the common emoji blocks.
# git grep exits 1 when it finds nothing, which under pipefail would abort
# this script precisely when the repository is clean. Hence "|| true".
emd=$( { git grep -In $'—' -- '*.md' '*.sh' '*.R' '*.py' '*.json' '*.tsv' 2>/dev/null || true; } | wc -l)
end=$( { git grep -In $'–' -- '*.md' '*.sh' '*.R' '*.py' '*.json' '*.tsv' 2>/dev/null || true; } | wc -l)
emo=$( { git grep -InP '[\x{1F300}-\x{1FAFF}\x{2600}-\x{27BF}\x{FE0F}]' \
        -- '*.md' '*.sh' '*.R' '*.py' 2>/dev/null || true; } | wc -l)
echo "  em dashes: ${emd}"
echo "  en dashes: ${end}"
echo "  emoji:     ${emo}"
[[ "${emd}" -eq 0 && "${emo}" -eq 0 ]] || { echo "  FAIL"; fail=1; }

note "4. banned phrases"
banned=("it is worth noting" "it is important to note" "in today's rapidly evolving"
        "plays a crucial role" "serves as a testament" "paving the way for"
        "in conclusion" "delve" "seamless" "underscores" "showcases"
        "not only" "leverage" "leverages" "leveraging")
hits=0
for b in "${banned[@]}"; do
    n=$( { git grep -Iiln -- "${b}" -- '*.md' '*.sh' '*.R' '*.py' 2>/dev/null || true; } | wc -l)
    [[ "${n}" -gt 0 ]] && { printf '  HIT  "%s" in %s file(s)\n' "${b}" "${n}"; hits=$((hits+1)); }
done
[[ "${hits}" -eq 0 ]] && echo "  none found"
[[ "${hits}" -eq 0 ]] || fail=1

note "5. nothing large tracked by git"
big=$(git ls-files -z | xargs -0 -r du -k 2>/dev/null | awk '$1 > 51200 {print $1" KB  "$2}')
if [[ -z "${big}" ]]; then echo "  no tracked file over 50 MB"; else echo "${big}"; fail=1; fi
echo "  largest tracked files:"
git ls-files -z | xargs -0 -r du -k 2>/dev/null | sort -rn | head -5 | awk '{printf "    %6d KB  %s\n", $1, $2}'

note "6. no data files tracked"
# grep -c exits 1 when it matches nothing, which is the expected case here.
tracked_data=$( { git ls-files || true; } | grep -Ec '\.(raw|RAW|mzML|mzML\.gz)$' || true)
echo "  RAW/mzML files tracked: ${tracked_data}"
[[ "${tracked_data}" -eq 0 ]] || fail=1

note "7. numbers and the file each comes from"
show() { [[ -f "$1" ]] && { printf '  --- %s\n' "$1"; sed 's/^/      /' "$1" | head -"${2:-12}"; } || printf '  (missing) %s\n' "$1"; }
show results/fasta_provenance.txt 14
show results/entrapment_build_stats.json 14
show results/maxquant_model_summary.tsv 9
show results/maxquant_filtering_waterfall.tsv 8
show results/maxquant_ribosomal_check.tsv 17
show results/maxquant_threshold_sensitivity.tsv 6
show results/entrapment_fdr_assessment.tsv 8
show results/engine_comparison.tsv 18
show results/conversion_sizes.tsv 8
show results/environment_versions.txt 20
show logs/stage_metrics.tsv 20

note "8. figures present and produced by a script"
for f in figures/*.png; do [[ -e "$f" ]] && printf '  %s\n' "$f"; done
echo "  drawn by scripts/figures.R via scripts/10_figures.sh"

printf '\n== result: %s ==\n' "$([[ "${fail}" -eq 0 ]] && echo PASS || echo "FAIL, see above")"
exit "${fail}"
