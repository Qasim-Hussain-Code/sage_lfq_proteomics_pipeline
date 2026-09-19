#!/usr/bin/env bash
# Build the two conda environments the pipeline needs and record exactly what
# versions landed. Nothing here pins a version taken from documentation: the
# script installs what the channels currently offer and then writes the
# resolved versions to results/environment_versions.txt. A search result is
# only reproducible if you can say which binary produced it.
#
# Two environments rather than one, because the R/Bioconductor stack and the
# mass spectrometry tools disagree about dependencies often enough that a
# single solve is slow and fragile.
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 01_install.sh [options]

Creates conda environments and records resolved versions.

Options:
  --ms-only       only the mass spectrometry environment (Sage, ThermoRawFileParser)
  --r-only        only the R/Bioconductor environment
  --force         recreate environments that already exist
  --dry-run       solve only, print what would be installed, change nothing
  --help          show this message
USAGE
}

MS_ONLY="false"; R_ONLY="false"; FORCE="false"; DRY_RUN="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ms-only) MS_ONLY="true"; shift ;;
        --r-only)  R_ONLY="true";  shift ;;
        --force)   FORCE="true";   shift ;;
        --dry-run) DRY_RUN="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "01_install: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${REPO_ROOT}/project.conf"
mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}"

# mamba resolves these environments in a fraction of the time conda takes and
# is a drop-in for the subcommands used here. Fall back if it is absent.
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; else SOLVER="conda"; fi
echo "01_install: using ${SOLVER}"

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

create_env() {
    local name="$1"; shift
    if env_exists "${name}" && [[ "${FORCE}" != "true" ]]; then
        echo "01_install: environment '${name}' already exists, skipping."
        return 0
    fi
    local args=(create --yes --name "${name}" "$@")
    [[ "${FORCE}" == "true" ]] && args=(create --yes --force --name "${name}" "$@")
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "01_install: dry run for '${name}'"
        "${SOLVER}" "${args[@]}" --dry-run
        return 0
    fi
    "${SOLVER}" "${args[@]}"
}

if [[ "${R_ONLY}" != "true" ]]; then
    # sage-proteomics carries the Sage binary. On the run that produced the
    # results in this repository the solver picked thermorawfileparser
    # 2.0.0.dev, which is the .NET 8 self-contained build, so no Mono was
    # pulled in and the environment came to 294 MB. I had expected the stable
    # 1.4.x line and a Mono dependency; pinning 1.4.x here would get that
    # instead. Both write the same gzipped mzML, so the rest of the pipeline
    # does not care which one is present.
    create_env "${CONDA_ENV_MS}" \
        --channel conda-forge --channel bioconda \
        sage-proteomics thermorawfileparser csvtk aria2 curl
fi

if [[ "${MS_ONLY}" != "true" ]]; then
    # QFeatures and msqrob2 come from bioconda's Bioconductor mirror. Building
    # these from source against a bare R takes a long time and needs a
    # compiler toolchain; the conda binaries do not.
    #
    # The Bioconductor packages are pinned and r-base is not. Asking for an
    # unpinned "r-base bioconductor-qfeatures bioconductor-msqrob2" instead
    # gave me a solve that quietly settled on R 4.3.3 and omitted QFeatures,
    # msqrob2 and MsCoreUtils altogether, leaving an environment that looked
    # installed and was missing everything that matters. Pinning the three
    # Bioconductor packages forces the solver to pick an R that can satisfy
    # them rather than an R that cannot. Versions are the bioconda latest as
    # of 2026-09-19; check for newer ones before reusing this.
    create_env "${CONDA_ENV_R}" \
        --channel conda-forge --channel bioconda \
        "bioconductor-qfeatures=${QFEATURES_VERSION:-1.20.0}" \
        "bioconductor-msqrob2=${MSQROB2_VERSION:-1.18.0}" \
        "bioconductor-mscoreutils=${MSCOREUTILS_VERSION:-1.22.1}" \
        bioconductor-limma bioconductor-biocparallel \
        r-ggplot2 r-data.table r-ggrepel r-patchwork
fi

[[ "${DRY_RUN}" == "true" ]] && { echo "01_install: dry run complete, nothing installed."; exit 0; }

# Record what actually landed. This file is the methods section for the
# software, and it is written from the installed binaries rather than from
# anything asserted above.
VERFILE="${RESULTS_DIR}/environment_versions.txt"
{
    echo "# Resolved tool versions"
    echo "# Written by scripts/01_install.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# Host: $(uname -srm)"
    echo
    if env_exists "${CONDA_ENV_MS}"; then
        echo "[${CONDA_ENV_MS}]"
        conda list --name "${CONDA_ENV_MS}" 2>/dev/null \
            | grep -iE '^(sage-proteomics|thermorawfileparser|mono|csvtk|aria2|dotnet)' || true
        echo
        echo "sage --version: $(conda run --name "${CONDA_ENV_MS}" sage --version 2>&1 | head -1)"
        echo "ThermoRawFileParser: $(conda run --name "${CONDA_ENV_MS}" ThermoRawFileParser --version 2>&1 | head -1)"
        echo
    fi
    if env_exists "${CONDA_ENV_R}"; then
        echo "[${CONDA_ENV_R}]"
        conda run --name "${CONDA_ENV_R}" Rscript -e \
            'cat(R.version.string, "\n"); for (p in c("QFeatures","msqrob2","MsCoreUtils","limma","BiocParallel","SummarizedExperiment","ggplot2")) cat(sprintf("%-22s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING")))' 2>/dev/null || true
    fi
} > "${VERFILE}"

echo "01_install: wrote ${VERFILE}"
cat "${VERFILE}"
