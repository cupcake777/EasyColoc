#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

usage() {
  cat <<'EOF'
Usage: ./easycoloc report-web RESULTS_DIR [--host HOST] [--port PORT] [--project-name NAME] [--refresh-data] [--open|--no-open]
EOF
}

host="127.0.0.1"
port="3000"
project_name=""
open_browser="false"
refresh_data="false"
results_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --host)
      host="${2:?missing value for --host}"
      shift 2
      ;;
    --port)
      port="${2:?missing value for --port}"
      shift 2
      ;;
    --project-name)
      project_name="${2:?missing value for --project-name}"
      shift 2
      ;;
    --refresh-data)
      refresh_data="true"
      shift
      ;;
    --open)
      open_browser="true"
      shift
      ;;
    --no-open)
      open_browser="false"
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    --*)
      echo "[REPORT-WEB] unexpected argument: $1" >&2
      exit 1
      ;;
    *)
      if [[ -z "${results_dir}" ]]; then
        results_dir="$1"
        shift
      else
        echo "[REPORT-WEB] unexpected positional argument: $1" >&2
        exit 1
      fi
      ;;
  esac
done

if [[ -z "${results_dir}" ]]; then
  echo "[REPORT-WEB] RESULTS_DIR is required" >&2
  usage >&2
  exit 1
fi

if [[ ! -d "${results_dir}" ]]; then
  echo "[REPORT-WEB] results directory not found: ${results_dir}" >&2
  exit 1
fi

if [[ ! -f "${repo_root}/web/dist/index.html" ]]; then
  echo "[REPORT-WEB] web/dist/index.html not found; run 'npm --prefix web run build' first" >&2
  exit 1
fi

data_file="${results_dir}/report_web/report-data.json"
if [[ "${refresh_data}" == "true" || ! -f "${data_file}" ]]; then
  build_args=(tools/build_report_web_data.R --results-dir "${results_dir}")
  if [[ -n "${project_name}" ]]; then
    build_args+=(--project-name "${project_name}")
  fi
  Rscript "${build_args[@]}"
fi

if [[ "${open_browser}" == "true" ]]; then
  echo "[REPORT-WEB] --open requested but browser auto-open is not enabled in this launcher"
fi

exec node tools/report_web_server.mjs \
  --host "${host}" \
  --port "${port}" \
  --data-file "${data_file}" \
  --dist-dir "${repo_root}/web/dist"
