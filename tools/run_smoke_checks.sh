#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "[CHECK] Parsing core scripts"
Rscript tests/check_parse.R

echo "[CHECK] Running plotting smoke test"
Rscript tests/smoke_test_plotting.R

echo "[CHECK] Running summary/report smoke test"
Rscript tests/smoke_test_summary_report.R

echo "[CHECK] Running runtime monitor smoke test"
Rscript tests/smoke_test_runtime_monitor.R

echo "[CHECK] Running doctor smoke test"
Rscript tests/smoke_test_doctor.R

echo "[CHECK] Running reference listing smoke test"
Rscript tests/smoke_test_refs.R

echo "[CHECK] Running reference bootstrap smoke test"
Rscript tests/smoke_test_bootstrap_refs.R

echo "[CHECK] Running toy demo bootstrap smoke test"
Rscript tests/smoke_test_bootstrap_demo.R

echo "[CHECK] Running 1KG setup smoke test"
Rscript tests/smoke_test_1kg_setup.R

echo "[CHECK] Running GTEx bootstrap smoke test"
Rscript tests/smoke_test_gtex_bootstrap.R

echo "[CHECK] Running LD cache smoke test"
Rscript tests/smoke_test_ld_cache.R

echo "[CHECK] Running CLI smoke test"
bash tests/smoke_test_cli.sh

echo "[CHECK] Running minimal synthetic demo"
Rscript examples/minimal/synthetic_plot_demo.R

echo "[CHECK] All smoke checks passed"
