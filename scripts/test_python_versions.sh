#!/usr/bin/env bash
# Testing script for maintainers:
# Test the current VAtools install ('pip install -e .`) across
# multiple Python versions using dedicated conda environments, creating
# any environment that doesn't already exist.
#
# Usage: scripts/test_python_versions.sh [python-version ...]
#   Defaults to: 3.10 3.11 3.12 3.13 3.14

set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_VERSIONS=("$@")
if [ ${#PYTHON_VERSIONS[@]} -eq 0 ]; then
    PYTHON_VERSIONS=("3.10" "3.11" "3.12" "3.13" "3.14")
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

LOG_DIR="$(mktemp -d)"
echo "Logs: ${LOG_DIR}"

declare -A STATUS
declare -A DETAIL

for ver in "${PYTHON_VERSIONS[@]}"; do
    env_name="vatools-py${ver}"
    install_log="${LOG_DIR}/install-${ver}.log"
    test_log="${LOG_DIR}/test-${ver}.log"

    echo ""
    echo "================================================================"
    echo "Python ${ver}  ->  conda env '${env_name}'"
    echo "================================================================"

    if conda env list | grep -qE "^${env_name}[[:space:]]"; then
        echo "Using existing conda env '${env_name}'"
    else
        echo "Creating conda env '${env_name}' (python=${ver})..."
        if ! conda create -y -n "${env_name}" "python=${ver}" > "${LOG_DIR}/create-${ver}.log" 2>&1; then
            STATUS[$ver]="ENV_CREATE_FAIL"
            DETAIL[$ver]="conda create failed; see ${LOG_DIR}/create-${ver}.log"
            echo "  -> FAILED to create environment"
            tail -n 15 "${LOG_DIR}/create-${ver}.log"
            continue
        fi
    fi

    conda activate "${env_name}"

    echo "--- pip install -e . ---"
    if pip install -e "${REPO_DIR}" > "${install_log}" 2>&1; then
        echo "--- python -m unittest discover ---"
        if (cd "${REPO_DIR}" && python -m unittest discover) > "${test_log}" 2>&1; then
            n_tests=$(grep -oE "Ran [0-9]+ test" "${test_log}" | grep -oE "[0-9]+" | head -1)
            STATUS[$ver]="PASS"
            DETAIL[$ver]="${n_tests:-?} tests passed"
            echo "  -> PASS (${DETAIL[$ver]})"
        else
            STATUS[$ver]="TEST_FAIL"
            DETAIL[$ver]="see ${test_log}"
            echo "  -> TEST FAILURES"
            tail -n 20 "${test_log}"
        fi
    else
        STATUS[$ver]="INSTALL_FAIL"
        DETAIL[$ver]="see ${install_log}"
        echo "  -> INSTALL FAILED"
        tail -n 20 "${install_log}"
    fi

    conda deactivate
done

echo ""
echo "================================================================"
echo "SUMMARY"
echo "================================================================"
overall_rc=0
for ver in "${PYTHON_VERSIONS[@]}"; do
    st="${STATUS[$ver]:-UNKNOWN}"
    printf "  Python %-6s : %-16s %s\n" "${ver}" "${st}" "${DETAIL[$ver]:-}"
    if [ "${st}" != "PASS" ]; then
        overall_rc=1
    fi
done
echo ""
echo "Full logs kept in: ${LOG_DIR}"
exit ${overall_rc}
