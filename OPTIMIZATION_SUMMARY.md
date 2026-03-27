# EasyColoc 优化实现总结

**日期:** 2026-03-27
**参考仓库:** ColocQuiaL, gwaslab

---

## 已完成的优化

### 1. APA 格式基因 ID 解析 (`parse_geneID`)

**位置:** `src/utils_format.R`

**功能:**
- 自动检测 GTEx APA 格式：`ENST001|GENE1|...`
- 支持标准 ENSEMBL ID (如 `ENSG00000123456.7`)
- 自动转换 ENSEMBL ID 到 HGNC 基因符号
- 返回结构化结果包含：
  - `geneSymbol`: 解析后的基因符号
  - `is_apa`: 是否为 APA 格式
  - `is_ensembl`: 是否为 ENSEMBL 格式
  - `transcript_id`: 转录本 ID (APA 格式)
  - `ensembl_id`: ENSEMBL ID (无版本号)

**使用示例:**
```r
# 单个基因 ID 解析
result <- parse_geneID("ENST001|ARL3|...")
print(result$geneSymbol)  # "ARL3"
print(result$is_apa)      # TRUE

# 批量处理
gene_ids <- c("ENSG00000123456", "ENST001|GENE1|...")
parsed <- apply_parse_geneID(gene_ids)
```

---

### 2. LiftOver 降级策略 (`run_liftover_fallback`)

**位置:** `src/utils_format.R`

**功能:**
- 当 gwaslab harmonization 失败时的 fallback 机制
- 迭代区域收缩策略（源自 ColocQuiaL）
- 参数:
  - `max_iterations`: 最大迭代次数 (默认 20)
  - `contraction_step`: 每次收缩碱基数 (默认 5000)
  - 最小区域阈值：10kb

**工作流程:**
1. 尝试完整区域的 LiftOver
2. 如果失败，从两端收缩区域
3. 重复直到成功或区域过小

**集成到 `run_gwaslab_harmonization`:**
```r
run_gwaslab_harmonization(
    sumstats_dt,
    ref_fasta,
    ...,
    liftover_chain = "/path/to/hg19ToHg38.over.chain.gz"  # 新增参数
)
```

---

### 3. 批量基因 ID 解析 (`apply_parse_geneID`)

**位置:** `src/utils_format.R`

**功能:**
- 批量应用 `parse_geneID` 到基因 ID 向量
- 返回统一格式的 data.frame

**使用示例:**
```r
qtl_genes <- c("ENSG00000123456", "ENSG00000123457")
parsed <- apply_parse_geneID(qtl_genes)
# 返回：original_id, geneSymbol, is_apa, is_ensembl, transcript_id, ensembl_id
```

---

## 现有功能的改进建议

### Gene Track 可视化
现有的 `genetrack()` 函数 (`src/utils_plot.R`) 已经实现完善：
- ✅ 支持 GTF 文件导入
- ✅ 支持 TxDb 注释回退
- ✅ 自动处理基因重叠（多层堆叠）
- ✅ 链特异性颜色编码
- ✅ ggrepel 标签防重叠

无需额外修改。

---

## 使用指南

### 在 run_coloc.r 中使用 APA 解析

在读取 QTL 数据后，处理基因 ID:
```r
# 假设 qtl_df 包含 gene_id 列
gene_info <- apply_parse_geneID(qtl_df$gene_id)
qtl_df$geneSymbol <- gene_info$geneSymbol
qtl_df$is_apa <- gene_info$is_apa
```

### 配置 LiftOver Fallback

在 `config/global.yaml` 中添加:
```yaml
harmonization_settings:
  env_name: "gwaslab"
  liftover_chain: "/path/to/hg19ToHg38.over.chain.gz"
```

---

## 测试计划

1. **APA 解析测试:**
   - GTEx v8 APA 数据格式
   - 标准 ENSEMBL ID
   - 无效 ID 处理

2. **LiftOver Fallback 测试:**
   - gwaslab 成功场景
   - gwaslab 失败 + fallback 成功
   - 两者均失败

---

## 参考

- **ColocQuiaL:** https://github.com/bychen9/ColocQuiaL
  - APA 格式解析逻辑
  - LiftOver 迭代收缩策略

- **gwaslab:** https://github.com/Cloufield/gwaslab
  - 数据标准化流程
  - 链特异性处理

---

## 后续优化建议

1. **增强 LD 参考匹配:**
   - 当 lead SNP 不在 LD 参考中时，自动选择替代 SNP
   - 已在 `LD_plot` 函数中部分实现

2. **Tabix 查询优化:**
   - 预筛选显著性位点减少查询量
   - 参考 ColocQuiaL 的 `create_tabix_*.sh` 脚本

3. **SuSiE 结果整合:**
   - 可信集 SNP 标注已在 `LD_plot` 中实现
   - 可考虑添加更多可视化选项
