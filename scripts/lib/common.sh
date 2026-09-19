#!/usr/bin/env bash
# Shared helpers. Sourced, never executed.
#
# The reason this file exists is the logging. The README reports measured peak
# resident memory and elapsed time per stage, and the only way that stays
# honest across a dozen scripts is if every stage records it the same way.

# measure_run <stage-label> <command...>
# Runs the command under GNU time, appends one row to logs/stage_metrics.tsv,
# and passes the command's exit status through unchanged. Falls back to plain
# execution with wall clock only when GNU time is unavailable, because busybox
# time and the bash builtin do not report maximum resident set size.
measure_run() {
    local label="$1"; shift
    local metrics="${LOGS_DIR:-logs}/stage_metrics.tsv"
    mkdir -p "$(dirname "${metrics}")"
    if [[ ! -s "${metrics}" ]]; then
        printf 'stage\tstarted_utc\telapsed_s\tpeak_rss_kb\texit_status\tcommand\n' > "${metrics}"
    fi
    local started tmp status elapsed rss
    started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    tmp="$(mktemp)"
    status=0
    if [[ -x /usr/bin/time ]]; then
        /usr/bin/time --format='%e %M' --output="${tmp}" "$@" || status=$?
        elapsed="$(awk '{print $1}' "${tmp}" | tail -1)"
        rss="$(awk '{print $2}' "${tmp}" | tail -1)"
    else
        local t0 t1
        t0="$(date +%s)"
        "$@" || status=$?
        t1="$(date +%s)"
        elapsed="$((t1 - t0))"
        rss="NA"
    fi
    rm -f "${tmp}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${label}" "${started}" "${elapsed:-NA}" "${rss:-NA}" "${status}" "$*" >> "${metrics}"
    return "${status}"
}

# free_gib <path> - free space on the filesystem holding path, in whole GiB.
free_gib() { df -BG --output=avail "$1" | tail -1 | tr -dc '0-9'; }

# die <message...>
die() { echo "error: $*" >&2; exit 1; }
