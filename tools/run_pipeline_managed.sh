#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

show_help() {
  cat <<'EOF'
Usage: tools/run_pipeline_managed.sh [--global PATH] [--gwas PATH] [--qtl PATH] [--output-dir PATH]

Starts EasyColoc with managed logging, runtime heartbeat, manifest build,
monitor snapshot and completion check.
EOF
}

global_cfg=""
gwas_cfg=""
qtl_cfg=""
output_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --global)
      global_cfg="${2:?missing value for --global}"
      shift 2
      ;;
    --gwas)
      gwas_cfg="${2:?missing value for --gwas}"
      shift 2
      ;;
    --qtl)
      qtl_cfg="${2:?missing value for --qtl}"
      shift 2
      ;;
    --output-dir)
      output_dir="${2:?missing value for --output-dir}"
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "[RUN] unknown argument: $1" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

if [[ -n "${global_cfg}" ]]; then
  export EASYCOLOC_GLOBAL_CONFIG="${global_cfg}"
fi
if [[ -n "${gwas_cfg}" ]]; then
  export EASYCOLOC_GWAS_CONFIG="${gwas_cfg}"
fi
if [[ -n "${qtl_cfg}" ]]; then
  export EASYCOLOC_QTL_CONFIG="${qtl_cfg}"
fi

if [[ -z "${output_dir}" ]]; then
  output_dir="$(Rscript -e 'source("src/utils_config.R"); cfg <- easycoloc_read_configs(); cat(normalizePath(cfg$global$output_dir, mustWork = FALSE))')"
else
  export EASYCOLOC_OUTPUT_DIR="${output_dir}"
fi

mkdir -p "${output_dir}"

timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="${output_dir}/run_easycoloc_full_${timestamp}.log"
run_label="managed_${timestamp}"

total_gwas="$(Rscript -e 'source("src/utils_config.R"); cfg <- easycoloc_read_configs(); cat(length(cfg$gwas$datasets))')"

echo "[RUN] log file: ${log_file}"
echo "[RUN] starting EasyColoc pipeline"

EASYCOLOC_LOG_FILE="${log_file}" \
EASYCOLOC_RUN_LABEL="${run_label}" \
EASYCOLOC_PARENT_PID="$$" \
Rscript run_coloc.R > "${log_file}" 2>&1 &
run_pid=$!

cleanup() {
  if kill -0 "${run_pid}" 2>/dev/null; then
    kill "${run_pid}" 2>/dev/null || true
  fi
}
trap cleanup INT TERM

render_progress() {
  local processed_gwas current_gwas abf_count rds_count plot_count filled bar_width
  processed_gwas="$(grep -c '^Processing GWAS Dataset:' "${log_file}" 2>/dev/null || true)"
  current_gwas="$(grep '^Processing GWAS Dataset:' "${log_file}" 2>/dev/null | tail -n 1 | sed 's/^Processing GWAS Dataset: //')"
  abf_count="$(find "${output_dir}/abf" -maxdepth 1 -type f -name '*_locus_results.csv' 2>/dev/null | wc -l | tr -d ' ')"
  rds_count="$(find "${output_dir}/rds" -maxdepth 1 -type f -name '*.rds' 2>/dev/null | wc -l | tr -d ' ')"
  plot_count="$(find "${output_dir}/plots" -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) 2>/dev/null | wc -l | tr -d ' ')"
  bar_width=24
  if [[ "${total_gwas}" -gt 0 ]]; then
    filled=$(( processed_gwas * bar_width / total_gwas ))
  else
    filled=0
  fi
  printf -v bar '%*s' "${filled}" ''
  bar="${bar// /#}"
  printf -v pad '%*s' $(( bar_width - filled )) ''
  pad="${pad// /-}"
  printf '\r[%-24s] GWAS %s/%s | current: %s | abf:%s rds:%s plots:%s' \
    "${bar}${pad}" "${processed_gwas}" "${total_gwas}" "${current_gwas:-starting}" "${abf_count}" "${rds_count}" "${plot_count}"
}

while kill -0 "${run_pid}" 2>/dev/null; do
  render_progress
  sleep 5
done

wait "${run_pid}"
exit_code=$?
printf '\n'

echo "[RUN] pipeline exited with code ${exit_code}"
echo "[RUN] building output manifest"
if Rscript tools/build_output_manifest.R "${output_dir}"; then
  echo "[RUN] manifest written"
else
  echo "[RUN] manifest build failed"
fi

echo "[RUN] collecting runtime monitor snapshot"
if Rscript tools/monitor_easycoloc.R "${output_dir}"; then
  echo "[RUN] monitor snapshot written"
else
  echo "[RUN] monitor snapshot failed"
fi

echo "[RUN] checking completion"
if Rscript tools/check_run_completion.R "${output_dir}"; then
  echo "[RUN] completion check passed"
else
  echo "[RUN] completion check failed"
fi

exit "${exit_code}"
