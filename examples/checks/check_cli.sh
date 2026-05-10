#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

bash -n easycoloc
bash -n tools/run_pipeline_managed.sh
bash -n tools/watch_output_dir.sh

tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

./easycoloc init "${tmp_dir}/project"

test -f "${tmp_dir}/project/config/global.yml"
test -f "${tmp_dir}/project/config/gwas.yml"
test -f "${tmp_dir}/project/config/qtl.yml"

./easycoloc refs \
  --global "${tmp_dir}/project/config/global.yml" \
  --gwas "${tmp_dir}/project/config/gwas.yml" \
  --qtl "${tmp_dir}/project/config/qtl.yml" >/dev/null

./easycoloc help | grep "report-web" >/dev/null
./easycoloc report-web --help >/dev/null

echo "[SMOKE] cli smoke test passed"
