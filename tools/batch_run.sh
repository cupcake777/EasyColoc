#!/usr/bin/env bash
# =============================================================================
# EasyColoc 批量处理 Wrapper 脚本
# =============================================================================
# 功能：一键运行多个 GWAS 数据集的 colocalization 分析
# 超越 ColocQuiaL 的单数据集限制
#
# 使用示例:
#   bash tools/batch_run.sh config/gwas_batch.yaml
#
# 批次文件格式 (config/gwas_batch.yaml):
#   datasets:
#     - id: "GWAS1"
#       file: "/path/to/gwas1.txt"
#       build: "hg19"
#       type: "cc"
#       prop: 0.5
#       columns: {...}
#     - id: "GWAS2"
#       file: "/path/to/gwas2.txt"
#       ...
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# 检查依赖
check_dependencies() {
    log_info "Checking dependencies..."

    local missing=()

    if ! command -v R &> /dev/null; then
        missing+=("R")
    fi

    if ! command -v python &> /dev/null; then
        missing+=("Python")
    fi

    if ! command -v plink &> /dev/null; then
        missing+=("PLINK")
    fi

    if ! command -v tabix &> /dev/null; then
        missing+=("Tabix")
    fi

    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Missing dependencies: ${missing[*]}"
        exit 1
    fi

    log_success "All dependencies found"
}

# 解析批次配置
parse_batch_config() {
    local config_file="$1"

    if [ ! -f "$config_file" ]; then
        log_error "Config file not found: $config_file"
        exit 1
    fi

    # 使用 Python 解析 YAML (更可靠)
    python3 -c "
import yaml
import json
import sys

try:
    with open('$config_file', 'r') as f:
        config = yaml.safe_load(f)

    datasets = config.get('datasets', [])
    if not datasets:
        print('ERROR: No datasets found in config', file=sys.stderr)
        sys.exit(1)

    # 输出为 JSON 供 bash 处理
    print(json.dumps(datasets))
except Exception as e:
    print(f'ERROR: {e}', file=sys.stderr)
    sys.exit(1)
"
}

# 运行单个 GWAS 数据集
run_single_gwas() {
    local gwas_id="$1"
    local gwas_file="$2"
    local gwas_build="$3"
    local gwas_type="$4"
    local gwas_prop="$5"
    local output_subdir="$6"

    log_info "========================================="
    log_info "Running GWAS: $gwas_id"
    log_info "========================================="

    # 创建临时配置
    local temp_gwas_config=$(mktemp)

    cat > "$temp_gwas_config" << EOF
datasets:
  - id: "$gwas_id"
    name: "$gwas_id"
    file: "$gwas_file"
    build: "$gwas_build"
    type: "$gwas_type"
    prop: ${gwas_prop:-0.5}
    columns:
      snp: "SNP"
      chrom: "CHR"
      pos: "POS"
      beta: "BETA"
      se: "SE"
      pval: "P"
      n: "N"
      maf: "EAF"
EOF

    # 复制全局配置
    cp "$PROJECT_DIR/config/global.yaml" "$output_subdir/config_global.yaml"
    cp "$PROJECT_DIR/config/qtl.yaml" "$output_subdir/config_qtl.yaml"
    cp "$temp_gwas_config" "$output_subdir/config_gwas.yaml"

    # 运行分析
    cd "$output_subdir"

    local log_file="$output_subdir/run_${gwas_id}.log"

    if Rscript "$PROJECT_DIR/run_coloc.r" >> "$log_file" 2>&1; then
        log_success "GWAS $gwas_id completed successfully"

        # 汇总结果
        if [ -f "$output_subdir/results/summary_statistics.txt" ]; then
            log_info "Summary for $gwas_id:"
            cat "$output_subdir/results/summary_statistics.txt"
        fi

        rm "$temp_gwas_config"
        return 0
    else
        log_error "GWAS $gwas_id failed. Check log: $log_file"
        rm "$temp_gwas_config"
        return 1
    fi
}

# 合并所有结果
merge_results() {
    local output_base="$1"
    local gwas_ids=("$@")

    log_info "========================================="
    log_info "Merging results from ${#gwas_ids[@]} GWAS datasets"
    log_info "========================================="

    # 使用 R 合并结果
    Rscript - << 'MERGE_SCRIPT'
library(data.table)
library(dplyr)

output_base <- Sys.getenv("OUTPUT_BASE", "/tmp")
gwas_ids <- strsplit(Sys.getenv("GWAS_IDS", ""), ",")[[1]]

all_results <- list()

for (gwas_id in gwas_ids) {
    result_file <- file.path(output_base, gwas_id, "results", "all_colocalization_results.csv")

    if (file.exists(result_file)) {
        dt <- fread(result_file)
        dt$GWAS_ID <- gwas_id
        all_results[[gwas_id]] <- dt
        message(sprintf("Loaded %d results from %s", nrow(dt), gwas_id))
    } else {
        warning(sprintf("No results found for %s", gwas_id))
    }
}

if (length(all_results) > 0) {
    merged <- rbindlist(all_results, fill = TRUE)

    # 按 PP4 排序
    merged <- merged[order(-PP4)]

    # 保存合并结果
    output_file <- file.path(output_base, "merged_colocalization_results.csv")
    fwrite(merged, output_file)
    message(sprintf("Saved merged results: %s (%d rows)", output_file, nrow(merged)))

    # 保存显著结果 (PP4 >= 0.7)
    sig_file <- file.path(output_base, "significant_colocalizations_PP4_0.7.csv")
    sig <- merged[PP4 >= 0.7]
    fwrite(sig, sig_file)
    message(sprintf("Saved significant results (PP4 >= 0.7): %s (%d rows)", sig_file, nrow(sig)))

    # 生成汇总统计
    summary_stats <- data.frame(
        Metric = c("Total Tests", "Significant (PP4>=0.7)", "Unique Loci",
                   "Unique Genes", "Mean PP4", "Max PP4"),
        Value = c(nrow(merged), nrow(sig),
                  length(unique(merged$Locus)),
                  length(unique(merged$Gene)),
                  round(mean(merged$PP4, na.rm = TRUE), 4),
                  round(max(merged$PP4, na.rm = TRUE), 4))
    )

    summary_file <- file.path(output_base, "batch_summary.csv")
    fwrite(summary_stats, summary_file)
    message(sprintf("Saved summary: %s", summary_file))

    # 生成汇总报告
    report_html <- file.path(output_base, "batch_report.html")

    sink(report_html)
    cat("<!DOCTYPE html>\n")
    cat("<html>\n<head>\n")
    cat("<title>EasyColoc Batch Analysis Report</title>\n")
    cat("<style>\n")
    cat("body { font-family: Arial, sans-serif; margin: 20px; }\n")
    cat("h1 { color: #2c3e50; }\n")
    cat("table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
    cat("th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n")
    cat("th { background-color: #3498db; color: white; }\n")
    cat("tr:nth-child(even) { background-color: #f2f2f2; }\n")
    cat(".highlight { background-color: #f1c40f; }\n")
    cat("</style>\n")
    cat("</head>\n<body>\n")
    cat(sprintf("<h1>EasyColoc Batch Analysis Report</h1>\n"))
    cat(sprintf("<p>Generated: %s</p>\n", Sys.time()))
    cat(sprintf("<p>Total GWAS datasets: %d</p>\n", length(gwas_ids)))
    cat("<h2>Summary Statistics</h2>\n")
    cat("<table>\n")
    cat("<tr><th>Metric</th><th>Value</th></tr>\n")
    for (i in 1:nrow(summary_stats)) {
        cat(sprintf("<tr><td>%s</td><td>%s</td></tr>\n",
                    summary_stats$Metric[i], summary_stats$Value[i]))
    }
    cat("</table>\n")

    if (nrow(sig) > 0) {
        cat("<h2>Top 20 Significant Colocalizations (PP4 >= 0.7)</h2>\n")
        cat("<table>\n")
        cat("<tr><th>GWAS</th><th>Locus</th><th>Gene</th><th>PP4</th><th>n_snps</th></tr>\n")
        top20 <- head(sig[order(-sig$PP4), ], 20)
        for (i in 1:nrow(top20)) {
            row_class <- if(top20$PP4[i] >= 0.8) 'class="highlight"' else ''
            cat(sprintf("<tr %s><td>%s</td><td>%s</td><td>%s</td><td>%.4f</td><td>%d</td></tr>\n",
                        row_class, top20$GWAS_ID[i], top20$Locus[i],
                        top20$Gene[i], top20$PP4[i], top20$n_snps[i]))
        }
        cat("</table>\n")
    }

    cat("</body>\n</html>\n")
    sink()

    message(sprintf("Saved HTML report: %s", report_html))
} else {
    stop("No results to merge")
}
MERGE_SCRIPT

    export OUTPUT_BASE="$output_base"
    export GWAS_IDS=$(IFS=,; echo "${gwas_ids[*]}")
}

# 主函数
main() {
    if [ $# -lt 1 ]; then
        echo "Usage: $0 <batch_config.yaml> [output_base]"
        echo ""
        echo "Example:"
        echo "  $0 config/gwas_batch.yaml"
        echo "  $0 config/gwas_batch.yaml /path/to/output"
        exit 1
    fi

    local batch_config="$1"
    local output_base="${2:-$PROJECT_DIR/batch_results}"

    log_info "========================================="
    log_info "EasyColoc Batch Processor"
    log_info "========================================="
    log_info "Config: $batch_config"
    log_info "Output: $output_base"

    check_dependencies

    # 创建输出目录
    mkdir -p "$output_base"

    # 解析批次配置
    log_info "Parsing batch configuration..."
    local datasets_json
    datasets_json=$(parse_batch_config "$batch_config")

    if [ $? -ne 0 ]; then
        log_error "Failed to parse batch config"
        exit 1
    fi

    # 获取数据集数量
    local num_datasets
    num_datasets=$(echo "$datasets_json" | python3 -c "import json,sys; print(len(json.load(sys.stdin)))")

    log_info "Found $num_datasets GWAS datasets to process"

    # 运行每个 GWAS 数据集
    local gwas_ids=()
    local failed_ids=()

    echo "$datasets_json" | python3 -c "
import json
import sys

datasets = json.load(sys.stdin)
for ds in datasets:
    print(f\"{ds.get('id', 'UNKNOWN')}|{ds.get('file', '')}|{ds.get('build', 'hg38')}|{ds.get('type', 'quant')}|{ds.get('prop', 0.5)}\")
" | while IFS='|' read -r gwas_id gwas_file gwas_build gwas_type gwas_prop; do

        gwas_ids+=("$gwas_id")

        local output_subdir="$output_base/$gwas_id"
        mkdir -p "$output_subdir"

        if ! run_single_gwas "$gwas_id" "$gwas_file" "$gwas_build" "$gwas_type" "$gwas_prop" "$output_subdir"; then
            failed_ids+=("$gwas_id")
        fi

    done

    # 合并结果
    if [ ${#gwas_ids[@]} -gt 0 ]; then
        merge_results "$output_base" "${gwas_ids[@]}"
    fi

    # 最终报告
    log_info "========================================="
    log_info "Batch Processing Complete"
    log_info "========================================="
    log_info "Output directory: $output_base"
    log_info "Processed: ${#gwas_ids[@]} datasets"

    if [ -f "$output_base/batch_report.html" ]; then
        log_success "View report: $output_base/batch_report.html"
    fi

    if [ ${#failed_ids[@]} -gt 0 ]; then
        log_warning "Failed datasets: ${failed_ids[*]}"
    fi
}

main "$@"
