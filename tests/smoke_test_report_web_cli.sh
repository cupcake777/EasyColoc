#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

tmp_dir="$(mktemp -d)"
trap 'if [[ -n "${server_pid:-}" ]]; then kill "${server_pid}" 2>/dev/null || true; fi; rm -rf "${tmp_dir}"' EXIT

mkdir -p "${tmp_dir}/results/abf" "${tmp_dir}/results/plots"
cat > "${tmp_dir}/results/abf/GWAS_A_rs100_locus_results.csv" <<'CSV'
GWAS_ID,QTL_ID,Locus,Phenotype,PP4,n_snps
GWAS_A,postnatal,rs100,GENE1,0.92,41
CSV

./easycoloc report-web "${tmp_dir}/results" --host 127.0.0.1 --port 43123 --no-open >"${tmp_dir}/report.log" 2>&1 &
server_pid=$!

for _ in $(seq 1 30); do
  if curl -fsS "http://127.0.0.1:43123/api/report-data" >/dev/null 2>&1; then
    break
  fi
  sleep 1
done

curl -fsS "http://127.0.0.1:43123/api/report-data" | grep '"project_name"' >/dev/null
grep '\[REPORT-WEB\] url:' "${tmp_dir}/report.log" >/dev/null

echo "[SMOKE] report-web cli smoke test passed"
