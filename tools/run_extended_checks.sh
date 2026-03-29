#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "[EXTENDED] Parsing core scripts"
Rscript tests/check_parse.R

echo "[EXTENDED] Running PLINK clump smoke test"
Rscript tests/smoke_test_plink_clump.R

echo "[EXTENDED] Running standard smoke checks"
./tools/run_smoke_checks.sh

echo "[EXTENDED] All extended checks passed"
