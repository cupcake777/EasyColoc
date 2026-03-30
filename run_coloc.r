#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
  library(parallel)
  library(tools)
  library(glue)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(jsonlite)
})
message("==============================================================")
message("EasyColoc v1.1 - Colocalization Analysis Pipeline")
message("==============================================================")

source("src/utils_config.R")

message("[INIT] Loading Utility Modules...")
utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

cfg_bundle <- easycoloc_read_configs()
cfg_global <- cfg_bundle$global
cfg_gwas   <- cfg_bundle$gwas
cfg_qtl    <- cfg_bundle$qtl
message(glue("[INIT] Configs: global={cfg_bundle$paths$global} gwas={cfg_bundle$paths$gwas} qtl={cfg_bundle$paths$qtl}"))

# =============================================================================
# Reproducibility: Set global random seed
# =============================================================================
# Seed controls: (1) Global R random number generation
#                (2) SuSiE-RSS iterative fine-mapping (uses offset seeds)
#                (3) Parallel mclapply (mc.set.seed = TRUE)
global_seed <- if(!is.null(cfg_global$random_seed)) {
    as.integer(cfg_global$random_seed)
} else {
    20240326  # Default seed if not specified
}
set.seed(global_seed)
message(glue("[INIT] Global random seed set: {global_seed}"))

# Temp directory for intermediate files (PLINK clump outputs, error logs, etc.)
temp_dir <- if (!is.null(cfg_global$temp_dir)) cfg_global$temp_dir else "./temp"

# Read colocalization analysis settings
coloc_pp4_thresh <- if(!is.null(cfg_global$coloc_settings$pp4_threshold)) {
    as.numeric(cfg_global$coloc_settings$pp4_threshold)
} else if(!is.null(cfg_global$sig_threshold)) {
    as.numeric(cfg_global$sig_threshold)
} else {
    0.8
}
susie_thresh <- if(!is.null(cfg_global$coloc_settings$susie_threshold)) {
    as.numeric(cfg_global$coloc_settings$susie_threshold)
} else {
    0.75
}
min_snps <- if(!is.null(cfg_global$coloc_settings$min_snps)) {
    as.integer(cfg_global$coloc_settings$min_snps)
} else {
    30
}
top_candidates <- if(!is.null(cfg_global$coloc_settings$top_candidates)) {
    as.integer(cfg_global$coloc_settings$top_candidates)
} else {
    100
}
maf_default <- if(!is.null(cfg_global$coloc_settings$maf_default)) {
    as.numeric(cfg_global$coloc_settings$maf_default)
} else {
    0.1
}
maf_na_replacement <- if(!is.null(cfg_global$coloc_settings$maf_na_replacement)) {
    as.numeric(cfg_global$coloc_settings$maf_na_replacement)
} else {
    0.05
}
maf_epsilon <- if(!is.null(cfg_global$coloc_settings$maf_epsilon)) {
    as.numeric(cfg_global$coloc_settings$maf_epsilon)
} else {
    1.0e-6
}
pvalue_floor <- if(!is.null(cfg_global$coloc_settings$pvalue_floor)) {
    as.numeric(cfg_global$coloc_settings$pvalue_floor)
} else {
    1.0e-300
}

# Read coloc priors from config
coloc_p1 <- if(!is.null(cfg_global$coloc_settings$p1)) as.numeric(cfg_global$coloc_settings$p1) else 1e-4
coloc_p2 <- if(!is.null(cfg_global$coloc_settings$p2)) as.numeric(cfg_global$coloc_settings$p2) else 1e-4
coloc_p12 <- if(!is.null(cfg_global$coloc_settings$p12)) as.numeric(cfg_global$coloc_settings$p12) else 5e-6

# Read harmonization settings

# Read harmonization settings
gwaslab_env <- if(!is.null(cfg_global$harmonization_settings$env_name)) {
    cfg_global$harmonization_settings$env_name
} else {
    "gwaslab"
}

# Read plot settings
plot_sig_threshold <- if(!is.null(cfg_global$plot_settings$significance_threshold)) {
    as.numeric(cfg_global$plot_settings$significance_threshold)
} else {
    5.0e-8
}
plot_window_bp <- if(!is.null(cfg_global$plot_settings$plot_window_bp)) {
    as.integer(cfg_global$plot_settings$plot_window_bp)
} else {
    200000
}
r2_breaks <- if(!is.null(cfg_global$plot_settings$r2_breaks)) {
    as.numeric(unlist(cfg_global$plot_settings$r2_breaks))
} else {
    c(0.2, 0.4, 0.6, 0.8, 1.0)
}
r2_colors <- if(!is.null(cfg_global$plot_settings$r2_colors)) {
    as.character(unlist(cfg_global$plot_settings$r2_colors))
} else {
    c("#313695", "#4575B4", "#74ADD1", "#FDB863", "#D73027")
}
lead_snp_color <- if(!is.null(cfg_global$plot_settings$lead_snp_color)) {
    cfg_global$plot_settings$lead_snp_color
} else {
    "#7F3C8D"
}
# Use clump p1 as significance threshold for plotting (or fall back to plot_settings)
clump_p1_for_plot <- if(!is.null(cfg_global$clump$p1)) as.numeric(cfg_global$clump$p1) else 5e-8
plot_sig_threshold <- if(!is.null(cfg_global$plot_settings$significance_threshold)) {
    as.numeric(cfg_global$plot_settings$significance_threshold)
} else {
    clump_p1_for_plot  # Default to clump threshold
}
# Generate significance label from threshold
significance_label <- if(!is.null(cfg_global$plot_settings$significance_label)) {
    cfg_global$plot_settings$significance_label
} else {
    # Auto-generate label from threshold (e.g., 5e-8 -> "5×10⁻⁸")
    thresh <- plot_sig_threshold
    if (thresh < 1e-6) {
        exp_val <- format(thresh, scientific = TRUE)
        # Convert "5e-08" to "5×10⁻⁸"
        exp_val <- gsub("e-0?", "×10⁻", exp_val)
        exp_val <- gsub("^-", "⁻", exp_val)
        paste0("GWAS ", exp_val)
    } else {
        # For larger thresholds, use simple format
        paste0("GWAS ", thresh)
    }
}

plot_width <- if(!is.null(cfg_global$plot_settings$plot_width)) {
    as.numeric(cfg_global$plot_settings$plot_width)
} else {
    10
}
plot_height <- if(!is.null(cfg_global$plot_settings$plot_height)) {
    as.numeric(cfg_global$plot_settings$plot_height)
} else {
    8
}
title_phenotype_field <- if(!is.null(cfg_global$plot_settings$title_phenotype_field)) {
    cfg_global$plot_settings$title_phenotype_field
} else {
    "gene"
}

# =============================================================================
# Analysis parameters (centralized for portability)
# =============================================================================
# flank_bp: Window size around lead SNP for QTL extraction
flank_bp <- if(!is.null(cfg_global$analysis$flank_bp)) {
    as.integer(cfg_global$analysis$flank_bp)
} else {
    500000  # Default 500kb window
}
prefilter_sig_pairs <- if (!is.null(cfg_global$analysis$prefilter_sig_pairs)) {
    isTRUE(cfg_global$analysis$prefilter_sig_pairs)
} else {
    TRUE
}
# min_snps_susie: Minimum SNPs required to run SuSiE fine-mapping
min_snps_susie <- if(!is.null(cfg_global$coloc_settings$min_snps_susie)) {
    as.integer(cfg_global$coloc_settings$min_snps_susie)
} else {
    10  # Default minimum for SuSiE
}

message(glue("[INIT] PP4 Threshold: {coloc_pp4_thresh}"))
message(glue("[INIT] SuSiE Threshold: {susie_thresh}"))
message(glue("[INIT] Min SNPs: {min_snps}"))
message(glue("[INIT] GWASLab Environment: {gwaslab_env}"))
message(glue("[INIT] QTL sigPairs prefilter: {prefilter_sig_pairs}"))

cfg_runtime <- if (!is.null(cfg_global$runtime)) cfg_global$runtime else list()
runtime_enabled <- if (!is.null(cfg_runtime$enabled)) isTRUE(cfg_runtime$enabled) else TRUE
resume_completed_tasks <- if (!is.null(cfg_runtime$resume_completed_tasks)) {
    isTRUE(cfg_runtime$resume_completed_tasks)
} else {
    TRUE
}
skip_existing_locus_results <- if (!is.null(cfg_runtime$skip_existing_locus_results)) {
    isTRUE(cfg_runtime$skip_existing_locus_results)
} else {
    FALSE
}
write_task_events <- if (!is.null(cfg_runtime$write_task_events)) {
    isTRUE(cfg_runtime$write_task_events)
} else {
    TRUE
}

base_out_dir <- normalizePath(cfg_global$output_dir, mustWork = FALSE)
dir_abf <- file.path(base_out_dir, "abf")
dir_susie <- file.path(base_out_dir, "susie")
dir_plots <- file.path(base_out_dir, "plots")
dir_rds <- file.path(base_out_dir, "rds")
for(d in c(base_out_dir, dir_abf, dir_susie, dir_plots, dir_rds)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

runtime_fingerprint <- compute_runtime_config_fingerprint(
    c(cfg_bundle$paths$global, cfg_bundle$paths$gwas, cfg_bundle$paths$qtl, "run_coloc.r")
)
initialize_runtime_tracker(
    output_dir = base_out_dir,
    enabled = runtime_enabled,
    config_fingerprint = runtime_fingerprint
)
register_active_run(
    log_file = Sys.getenv("EASYCOLOC_LOG_FILE", unset = NA_character_),
    run_label = Sys.getenv("EASYCOLOC_RUN_LABEL", unset = NA_character_),
    parent_pid = suppressWarnings(as.integer(Sys.getenv("EASYCOLOC_PARENT_PID", unset = NA_character_))),
    command = "Rscript run_coloc.r"
)
on.exit(clear_active_run(), add = TRUE)
append_runtime_event(
    stage = "init",
    message_text = "EasyColoc pipeline initialized",
    extras = list(output_dir = base_out_dir, random_seed = global_seed)
)
write_runtime_heartbeat(
    stage = "init",
    message_text = "Pipeline initialized",
    counters = list(total_gwas = length(cfg_gwas$datasets))
)

message(glue("[INIT] Output: {base_out_dir}"))
hash_table <- TRUE
if (!is.null(cfg_global$hash_table_dir) && dir.exists(cfg_global$hash_table_dir)) {
    hash_table <- initialize_hash_system(cfg_global$hash_table_dir)
}

identify_loci <- function(sumstats_dt, p_col, snp_col, chrom_col, pos_col, plink_bfile = NULL, plink_bin = "plink",
                           clump_p1 = 5e-8, clump_p2 = 5e-8, clump_kb = 1000, clump_r2 = 0.1, dataset_id = NULL, keep_file = NULL,
                           ea_col = NULL, nea_col = NULL) {
    message(glue("[LOCUS] Identifying loci (P < {clump_p1})..."))
    if (!all(c(p_col, snp_col, chrom_col, pos_col) %in% names(sumstats_dt))) {
        stop("Columns not found in sumstats_dt")
    }
    sig_dt <- as.data.table(copy(sumstats_dt[sumstats_dt[[p_col]] <= clump_p1, ]))
    if (nrow(sig_dt) == 0) {
        warning("No significant SNPs.")
        return(NULL)
    }

    if (is.null(plink_bfile)) {
        stop("[LOCUS] PLINK bfile is required for clumping")
    }
    if (is.null(ea_col) || is.null(nea_col) || !all(c(ea_col, nea_col) %in% names(sig_dt))) {
        stop("[LOCUS] EA/NEA columns are required for PLINK ID translation")
    }

    normalize_chr_for_plink <- function(x) {
        gsub("^chr", "", as.character(x), ignore.case = TRUE)
    }

    message(glue("        Using PLINK clumping with: {basename(plink_bfile)}"))
    if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
    ds_id <- if (is.null(dataset_id)) "unknown" else dataset_id
    temp_assoc <- file.path(temp_dir, glue("clump_input_{ds_id}.qassoc"))
    temp_prefix <- file.path(temp_dir, glue("clump_output_{ds_id}"))

    ref_panel <- load_plink_bim_index(plink_bfile)
    if (is.null(ref_panel) || nrow(ref_panel) == 0) {
        stop(glue("[LOCUS] Failed to load PLINK BIM index for {plink_bfile}"))
    }

    sig_dt[, CHR_CLUMP := normalize_chr_for_plink(get(chrom_col))]
    sig_dt[, POS_CLUMP := suppressWarnings(as.numeric(get(pos_col)))]
    sig_dt[, EA_CLUMP := as.character(get(ea_col))]
    sig_dt[, NEA_CLUMP := as.character(get(nea_col))]

    matched_dt <- merge(
        sig_dt,
        ref_panel,
        by.x = c("CHR_CLUMP", "POS_CLUMP"),
        by.y = c("CHR", "POS"),
        all.x = FALSE,
        all.y = FALSE
    )

    matched_dt <- matched_dt[
        (EA_CLUMP == V5 & NEA_CLUMP == V6) |
        (EA_CLUMP == V6 & NEA_CLUMP == V5)
    ]

    if (nrow(matched_dt) == 0) {
        stop(glue("[LOCUS] No variants matched the PLINK reference after CHR/POS + dual-direction allele matching for {ds_id}"))
    }

    message(glue("        Matched {nrow(matched_dt)} GWAS rows to PLINK SNP IDs"))

    setorderv(matched_dt, cols = p_col, order = 1L, na.last = TRUE)
    matched_dt <- matched_dt[!duplicated(plink_snp_id)]

    clump_write <- matched_dt[, .(SNP = plink_snp_id, P = get(p_col))]
    clump_write <- clump_write[!is.na(P) & P > 0]

    if (nrow(clump_write) == 0) {
        stop(glue("[LOCUS] No valid PLINK clump input rows remain for {ds_id} after ID translation"))
    }

    fwrite(clump_write, temp_assoc, sep = "\t", quote = FALSE, col.names = TRUE)
    message(glue("        Wrote {nrow(clump_write)} translated SNPs to {temp_assoc}"))

    keep_cmd <- if (!is.null(keep_file) && file.exists(keep_file)) {
        glue("--keep {keep_file}")
    } else {
        ""
    }

    cmd <- glue(
        "{plink_bin} --bfile {plink_bfile} --clump {temp_assoc} {keep_cmd} ",
        "--clump-p1 {clump_p1} --clump-p2 {clump_p2} ",
        "--clump-r2 {clump_r2} --clump-kb {clump_kb} --out {temp_prefix}"
    )

    message("        Running PLINK clumping...")
    sys_result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    log_file <- paste0(temp_prefix, ".log")
    if (file.exists(log_file)) {
        log_content <- readLines(log_file, warn = FALSE)
        clump_msg <- grep("clumps formed|No significant --clump results", log_content, value = TRUE)
        if (length(clump_msg) > 0) {
            message(glue("        {trimws(clump_msg[length(clump_msg)])}"))
        }
    }

    if (sys_result != 0) {
        stop(glue("[LOCUS] PLINK execution failed for {ds_id} (exit code: {sys_result})"))
    }

    clump_file <- paste0(temp_prefix, ".clumped")
    if (!file.exists(clump_file)) {
        stop(glue("[LOCUS] PLINK did not produce a .clumped file for {ds_id}"))
    }

    clump_res <- fread(clump_file)
    if (nrow(clump_res) == 0) {
        warning(glue("[LOCUS] PLINK produced 0 clumps for {ds_id} after rsID translation"))
        return(list())
    }

    message(glue("[LOCUS] Found {nrow(clump_res)} clumps via PLINK."))

    loci <- vector("list", nrow(clump_res))
    for (i in seq_len(nrow(clump_res))) {
        lead_plink_id <- clump_res$SNP[i]
        row <- matched_dt[plink_snp_id == lead_plink_id][1]
        if (nrow(row) == 0) next

        loci[[i]] <- list(
            chrom = as.character(row[[chrom_col]]),
            pos = as.numeric(row[[pos_col]]),
            snp = row[[snp_col]],
            p = as.numeric(row[[p_col]]),
            identification_method = "PLINK_clump"
        )
    }

    Filter(Negate(is.null), loci)
}


run_pipeline <- function() {
    qtl_meta <- fread(cfg_qtl$qtl_info$file)
    qtl_meta <- easycoloc_resolve_qtl_meta_paths(qtl_meta, cfg_qtl$qtl_info$file, cfg_qtl)
    message(glue("[DATA] Loaded {nrow(qtl_meta)} QTL datasets."))

    if (!is.null(cfg_global$n_cores)) {
        n_cores <- as.integer(cfg_global$n_cores)
    } else if (isTRUE(cfg_global$use_parallel)) {
        n_cores <- parallel::detectCores() - 1
    } else {
        n_cores <- 1
    }
    n_cores <- max(1L, n_cores)
    qtl_parallel_min <- if (!is.null(cfg_global$runtime$qtl_parallel_min)) {
        max(1L, as.integer(cfg_global$runtime$qtl_parallel_min))
    } else {
        4L
    }
    # Forking across QTL tasks interacts poorly with SuSiE, LD extraction and
    # graphics devices; keep per-locus task execution serial unless explicitly
    # enabled in runtime settings.
    parallel_qtl_tasks <- if (!is.null(cfg_global$runtime$parallel_qtl_tasks)) {
        isTRUE(cfg_global$runtime$parallel_qtl_tasks)
    } else {
        FALSE
    }
    total_gwas <- length(cfg_gwas$datasets)
    message(glue("[INIT] Using {n_cores} CPU cores"))
    message(glue("[INIT] Per-locus QTL parallelism: {if (parallel_qtl_tasks) 'enabled' else 'disabled'}"))
    write_runtime_heartbeat(
        stage = "pipeline_start",
        message_text = "Main GWAS loop started",
        counters = list(total_gwas = total_gwas, total_qtl = nrow(qtl_meta), n_cores = n_cores)
    )

    for (gwas_idx in seq_along(cfg_gwas$datasets)) {
        gwas_cfg <- cfg_gwas$datasets[[gwas_idx]]
        separator <- strrep("=", 50)
        message(glue("\n{separator}\nProcessing GWAS Dataset: {gwas_cfg$name}\n{separator}"))
        append_runtime_event(
            stage = "gwas_start",
            message_text = "Processing GWAS dataset",
            gwas_id = gwas_cfg$id,
            extras = list(gwas_index = gwas_idx, total_gwas = total_gwas)
        )
        write_runtime_heartbeat(
            stage = "gwas_start",
            message_text = gwas_cfg$id,
            counters = list(gwas_index = gwas_idx, total_gwas = total_gwas)
        )

        tryCatch({
            if (!file.exists(gwas_cfg$file)) {
                warning(glue("File missing: {gwas_cfg$file}"))
                append_runtime_event(
                    level = "WARN",
                    stage = "gwas_missing",
                    message_text = "GWAS input file missing",
                    gwas_id = gwas_cfg$id,
                    extras = list(file = gwas_cfg$file)
                )
            } else {
                gwas_select_cols <- unique(unname(unlist(gwas_cfg$columns, use.names = FALSE)))
                gwas_select_cols <- gwas_select_cols[!is.na(gwas_select_cols) & nzchar(gwas_select_cols)]
                gwas_raw <- fread(gwas_cfg$file, select = gwas_select_cols, showProgress = FALSE)
                gwas_std <- format_sumstats(
                    gwas_raw,
                    type = "gwas",
                    col_map = as.list(gwas_cfg$columns),
                    case_control = (gwas_cfg$type == "cc")
                )

                ref_fasta <- if (gwas_cfg$build == "hg19") cfg_global$ref_genome_hg19 else cfg_global$ref_genome_hg38
                pop <- if (is.null(gwas_cfg$pop)) "EUR" else gwas_cfg$pop
                ref_vcf_1kg <- file.path(
                    cfg_global[["1kg_af"]],
                    glue("1KG_hg{gsub('hg', '', gwas_cfg$build)}_{pop}_AF.tsv.gz")
                )
                ref_dbsnp <- if (gwas_cfg$build == "hg19") cfg_global[["dbsnp_hg19"]] else cfg_global[["dbsnp_hg38"]]
                ref_alt_freq <- "AF"
                qtl_build <- if (is.null(cfg_qtl$qtl_info$build)) "hg38" else cfg_qtl$qtl_info$build

                gwas_harm <- run_gwaslab_harmonization(
                    gwas_std,
                    ref_fasta = ref_fasta,
                    ref_vcf = ref_vcf_1kg,
                    ref_dbsnp = ref_dbsnp,
                    ref_alt_freq = ref_alt_freq,
                    source_build = if (is.null(gwas_cfg$build)) "19" else gsub("hg", "", gwas_cfg$build),
                    target_build = gsub("hg", "", qtl_build),
                    n_threads = cfg_global$n_cores,
                    save_dir = cfg_global$harmonize_dir,
                    dataset_id = gwas_cfg$id,
                    env_name = gwaslab_env
                )

                if (!"N" %in% names(gwas_harm) && !is.null(gwas_cfg$sample_size_n)) {
                    gwas_harm[, N := as.numeric(gwas_cfg$sample_size_n)]
                }

                plink_build <- if (is.null(cfg_qtl$qtl_info$build)) "hg38" else as.character(cfg_qtl$qtl_info$build)
                plink_ref_hg38 <- if (identical(plink_build, "hg19")) cfg_global$plink_hg19 else cfg_global$plink_hg38
                if (is.null(plink_ref_hg38) || !nzchar(plink_ref_hg38)) {
                    stop(glue("PLINK reference not configured for QTL build: {plink_build}"))
                }
                rsid_col_gwas <- if ("rsID" %in% names(gwas_harm)) "rsID" else if ("SNPID" %in% names(gwas_harm)) "SNPID"
                p_col_locus <- if ("P" %in% names(gwas_harm)) "P"
                chr_col_locus <- if ("CHR" %in% names(gwas_harm)) "CHR"
                pos_col_locus <- if ("POS" %in% names(gwas_harm)) "POS"
                if (is.null(p_col_locus) || is.null(rsid_col_gwas) || is.null(chr_col_locus) || is.null(pos_col_locus)) {
                    stop("Missing required columns for locus identification")
                }

                clump_p1 <- if (!is.null(gwas_cfg$clump_p1)) {
                    as.numeric(gwas_cfg$clump_p1)
                } else if (!is.null(gwas_cfg$clump_thres)) {
                    as.numeric(gwas_cfg$clump_thres)
                } else if (!is.null(cfg_global$clump$p1)) {
                    as.numeric(cfg_global$clump$p1)
                } else {
                    5e-8
                }

                clump_p2 <- if (!is.null(gwas_cfg$clump_p2)) {
                    as.numeric(gwas_cfg$clump_p2)
                } else if (!is.null(cfg_global$clump$p2)) {
                    as.numeric(cfg_global$clump$p2)
                } else {
                    clump_p1
                }

                clump_kb <- if (!is.null(gwas_cfg$clump_kb)) {
                    as.numeric(gwas_cfg$clump_kb)
                } else if (!is.null(cfg_global$clump$kb)) {
                    as.numeric(cfg_global$clump$kb)
                } else {
                    1000
                }

                clump_r2 <- if (!is.null(gwas_cfg$clump_r2)) {
                    as.numeric(gwas_cfg$clump_r2)
                } else if (!is.null(cfg_global$clump$r2)) {
                    as.numeric(cfg_global$clump$r2)
                } else {
                    0.1
                }

                message(glue("[LOCUS] Clump params: p1={clump_p1}, p2={clump_p2}, kb={clump_kb}, r2={clump_r2}"))

                keep_file <- NULL
                if (!is.null(cfg_global$plink_keep) && nzchar(cfg_global$plink_keep)) {
                    if (file.exists(cfg_global$plink_keep)) {
                        keep_file <- cfg_global$plink_keep
                    } else {
                        keep_base_dir <- dirname(cfg_global$plink_keep)
                        candidate_keep <- file.path(keep_base_dir, glue("{pop}.sample"))
                        if (file.exists(candidate_keep)) keep_file <- candidate_keep
                    }
                }
                if (is.null(keep_file)) {
                    message(glue("        WARNING: Keep file not found for {pop}, using all samples"))
                } else {
                    message(glue("        Using population-specific LD: {pop}"))
                }

                ea_col_locus <- if ("EA" %in% names(gwas_harm)) "EA" else NULL
                nea_col_locus <- if ("NEA" %in% names(gwas_harm)) "NEA" else NULL

                loci_list <- identify_loci(
                    gwas_harm,
                    p_col = p_col_locus,
                    snp_col = rsid_col_gwas,
                    chrom_col = chr_col_locus,
                    pos_col = pos_col_locus,
                    plink_bfile = plink_ref_hg38,
                    clump_p1 = clump_p1,
                    clump_p2 = clump_p2,
                    clump_kb = clump_kb,
                    clump_r2 = clump_r2,
                    dataset_id = gwas_cfg$id,
                    keep_file = keep_file,
                    ea_col = ea_col_locus,
                    nea_col = nea_col_locus
                )

                if (is.null(loci_list) || length(loci_list) == 0) {
                    message("No significant loci.")
                    append_runtime_event(
                        stage = "gwas_no_loci",
                        message_text = "No significant loci found",
                        gwas_id = gwas_cfg$id
                    )
                } else {
                    for (locus_idx in seq_along(loci_list)) {
                        locus <- loci_list[[locus_idx]]
                        message(glue("\n>>> Locus: {locus$snp} (hg38 chr{locus$chrom}:{locus$pos})"))
                        query_chrom <- locus$chrom
                query_pos <- locus$pos
                locus_id <- glue("chr{query_chrom}:{query_pos}:{locus$snp}")
                qtl_start <- max(1, query_pos - flank_bp)
                qtl_end <- query_pos + flank_bp
                message(glue("[LOCUS] Extraction window: ±{flank_bp/1000}kb ({qtl_start}-{qtl_end})"))

                append_runtime_event(
                    stage = "locus_start",
                    message_text = "Processing locus",
                    gwas_id = gwas_cfg$id,
                    locus_id = locus_id,
                    extras = list(locus_index = locus_idx, total_loci = length(loci_list))
                )
                write_runtime_heartbeat(
                    stage = "locus_start",
                    message_text = locus_id,
                    counters = list(
                        gwas_index = gwas_idx,
                        total_gwas = total_gwas,
                        locus_index = locus_idx,
                        total_loci = length(loci_list)
                    )
                )

                gwas_locus <- as.data.table(gwas_harm[CHR == query_chrom & POS >= qtl_start & POS <= qtl_end, ])
                if (nrow(gwas_locus) == 0) next
                if (!is.null(rsid_col_gwas) && rsid_col_gwas %in% names(gwas_locus)) {
                    gwas_locus[, rsid := get(rsid_col_gwas)]
                }
                locus_ld_res <- NULL
                if ("rsid" %in% names(gwas_locus)) {
                    locus_rsids <- unique(as.character(gwas_locus$rsid))
                    locus_rsids <- locus_rsids[!is.na(locus_rsids) & grepl("^rs", locus_rsids)]
                    if (length(locus_rsids) >= min_snps_susie) {
                        locus_ld_res <- get_ld_matrix(
                            locus_rsids,
                            plink_ref_hg38,
                            "plink",
                            keep_file = keep_file
                        )
                    }
                }

                recomb_data <- load_recomb_map(query_chrom, qtl_start, qtl_end, cfg_global$recom)
                results_fname <- sanitize_filename(glue("{gwas_cfg$id}_{locus$snp}_locus_results.csv"))
                results_path <- file.path(dir_abf, results_fname)

                if (resume_completed_tasks &&
                    skip_existing_locus_results &&
                    file.exists(results_path) &&
                    !is.na(file.info(results_path)$size) &&
                    file.info(results_path)$size > 0) {
                    message(glue("[RESUME] Existing locus result found, skipping: {basename(results_path)}"))
                    append_runtime_event(
                        stage = "locus_skipped_existing",
                        message_text = "Skipped locus because result file already exists",
                        gwas_id = gwas_cfg$id,
                        locus_id = locus_id,
                        extras = list(result_file = results_path)
                    )
                    next
                }

                process_qtl_wrapper <- function(row_idx) {
                    meta_row <- qtl_meta[row_idx, ]
                    qtl_id <- as.character(meta_row[[cfg_qtl$qtl_info$columns$id]])
                    qtl_file <- meta_row[[cfg_qtl$qtl_info$columns$all_filename]]
                    qtl_n <- as.numeric(meta_row[[cfg_qtl$qtl_info$columns$sample_size]])
                    candidate_phenos <- get_qtl_candidate_phenotypes(
                        meta_row = meta_row,
                        cfg_qtl = cfg_qtl,
                        chrom = query_chrom,
                        start = qtl_start,
                        end = qtl_end,
                        prefilter_sig_pairs = prefilter_sig_pairs,
                        verbose = FALSE
                    )

                    if (length(candidate_phenos) == 0 && isTRUE(prefilter_sig_pairs)) {
                        return(list(records = list(), rows = list()))
                    }

                    qtl_raw <- tryCatch(
                        query_tabix_region(qtl_file, query_chrom, qtl_start, qtl_end),
                        error = function(e) NULL
                    )
                    if (is.null(qtl_raw) || nrow(qtl_raw) == 0) {
                        return(list(records = list(), rows = list()))
                    }

                    if (ncol(qtl_raw) == length(cfg_qtl$QTL_all_header)) {
                        colnames(qtl_raw) <- cfg_qtl$QTL_all_header
                    }

                    pheno_col <- cfg_qtl$QTL_cols$phenotype
                    if (!is.null(candidate_phenos) &&
                        length(candidate_phenos) > 0 &&
                        !is.null(pheno_col) &&
                        pheno_col %in% names(qtl_raw)) {
                        qtl_raw <- qtl_raw[qtl_raw[[pheno_col]] %in% candidate_phenos, ]
                        if (nrow(qtl_raw) == 0) {
                            return(list(records = list(), rows = list()))
                        }
                    }
                    qtl_std_all <- format_sumstats(qtl_raw, type = "qtl", col_map = as.list(cfg_qtl$QTL_cols))
                    if ("variant_id" %in% names(qtl_std_all)) qtl_std_all[, rsid := variant_id]

                    qtl_groups <- if (is.null(pheno_col) || !pheno_col %in% names(qtl_std_all)) {
                        list(Combined = seq_len(nrow(qtl_std_all)))
                    } else {
                        split(seq_len(nrow(qtl_std_all)), qtl_std_all[[pheno_col]])
                    }
                    unique_phenos <- names(qtl_groups)
                    unique_phenos <- unique_phenos[!is.na(unique_phenos) & nzchar(unique_phenos)]

                    records <- list()
                    rows <- list()

                    for (pheno in unique_phenos) {
                        task_id <- build_task_id(gwas_cfg$id, locus_id, qtl_id, pheno)

                        if (resume_completed_tasks && runtime_task_completed(task_id, runtime_fingerprint)) {
                            prev_task <- runtime_get_task_record(task_id)
                            if (!is.null(prev_task)) {
                                rows[[length(rows) + 1]] <- data.table(
                                    GWAS_ID = prev_task$gwas_id,
                                    QTL_ID = prev_task$qtl_id,
                                    Locus = locus$snp,
                                    Phenotype = prev_task$phenotype,
                                    PP4 = suppressWarnings(as.numeric(prev_task$pp4)),
                                    n_snps = suppressWarnings(as.numeric(prev_task$n_snps)),
                                    identification_method = if ("identification_method" %in% names(prev_task)) {
                                        prev_task$identification_method
                                    } else {
                                        locus$identification_method
                                    }
                                )
                            }
                            records[[length(records) + 1]] <- list(
                                task_id = task_id,
                                status = "skipped",
                                gwas_id = gwas_cfg$id,
                                locus_id = locus_id,
                                qtl_id = qtl_id,
                                phenotype = pheno,
                                identification_method = locus$identification_method,
                                pp4 = if (!is.null(prev_task)) prev_task$pp4 else NA_real_,
                                n_snps = if (!is.null(prev_task)) prev_task$n_snps else NA_real_,
                                result_file = results_path,
                                plot_file = if (!is.null(prev_task) && "plot_file" %in% names(prev_task)) prev_task$plot_file else NA_character_,
                                rds_file = if (!is.null(prev_task) && "rds_file" %in% names(prev_task)) prev_task$rds_file else NA_character_,
                                error_message = NA_character_
                            )
                            next
                        }

                        task_result <- tryCatch({
                            qtl_idx <- qtl_groups[[pheno]]
                            qtl_std <- qtl_std_all[qtl_idx, ]

                            input_data <- prep_coloc_input_file(
                                gwas_locus, qtl_std,
                                list(pval = "P", beta = "BETA", se = "SE", n = "N"),
                                list(pval = "P", beta = "BETA", se = "SE"),
                                use_hash_table = hash_table,
                                min_snps = min_snps,
                                pvalue_floor = pvalue_floor
                            )

                            if (nrow(input_data) < min_snps) {
                                return(list(
                                    record = list(
                                        task_id = task_id,
                                        status = "filtered",
                                        gwas_id = gwas_cfg$id,
                                        locus_id = locus_id,
                                        qtl_id = qtl_id,
                                        phenotype = pheno,
                                        identification_method = locus$identification_method,
                                        pp4 = NA_real_,
                                        n_snps = nrow(input_data),
                                        result_file = results_path,
                                        plot_file = NA_character_,
                                        rds_file = NA_character_,
                                        error_message = NA_character_
                                    ),
                                    row = NULL
                                ))
                            }

                            res <- get_coloc_results(
                                input_data, gwas_cfg$type, gwas_cfg$prop, gwas_cfg$sample_size_n, qtl_n,
                                plink_bfile = plink_ref_hg38,
                                plink_bin = "plink",
                                use_susie = TRUE,
                                susie_threshold = susie_thresh,
                                maf_default = maf_default,
                                maf_na_replacement = maf_na_replacement,
                                maf_epsilon = maf_epsilon,
                                keep_file = keep_file,
                                ld_res = locus_ld_res,
                                p1 = coloc_p1,
                                p2 = coloc_p2,
                                p12 = coloc_p12
                            )

                            pp4 <- as.numeric(res$summary["PP.H4.abf"])
                            nsnps <- as.numeric(res$summary["nsnps"])
                            gene_sym_for_filename <- resolve_gene_label(pheno)
                            plot_path <- NA_character_
                            rds_path <- NA_character_

                            if (!is.na(pp4) && pp4 > coloc_pp4_thresh) {
                                if (!is.null(res$susie_result) &&
                                    !is.null(res$susie_result$summary) &&
                                    nrow(res$susie_result$summary) > 0) {
                                    fname <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_susie.csv"))
                                    fwrite(as.data.frame(res$susie_result$summary), file.path(dir_susie, fname))
                                }
                            }

                            if (!is.na(pp4) && pp4 > coloc_pp4_thresh) {
                                credible_set_snps <- NULL
                                if (!is.null(res$results) && "SNP.PP.H4" %in% names(res$results)) {
                                    snp_pp <- res$results[order(-res$results$SNP.PP.H4), c("snp", "SNP.PP.H4")]
                                    snp_pp <- snp_pp[snp_pp$SNP.PP.H4 > 0, ]
                                    if (nrow(snp_pp) > 0) {
                                        cumsum_pp <- cumsum(snp_pp$SNP.PP.H4)
                                        n_in_set <- min(which(cumsum_pp >= 0.95), nrow(snp_pp))
                                        credible_set_snps <- snp_pp[1:n_in_set, ]
                                        message(glue("  -> 95% credible set: {nrow(credible_set_snps)} SNP(s) - top: {credible_set_snps$snp[1]} (PP.H4={round(credible_set_snps$SNP.PP.H4[1], 4)})"))
                                    }
                                }

                                tryCatch({
                                    assign("lead_SNP", locus$snp, envir = .GlobalEnv)
                                    assign("geneSymbol", gene_sym_for_filename, envir = .GlobalEnv)
                                    assign("plink_bfile", plink_ref_hg38, envir = .GlobalEnv)
                                    assign("trait", gwas_cfg$id, envir = .GlobalEnv)

                                    p <- plot_qtl_association(
                                        qtl_all_chrom = "CHR.qtl",
                                        qtl_all_pvalue = "P.qtl",
                                        leadSNP_DF = input_data,
                                        ld_df = NULL,
                                        gtf_path = cfg_global$gene_anno,
                                        region_recomb = recomb_data,
                                        recomb_path = cfg_global$recom,
                                        show_lead_line = FALSE,
                                        qtl_type = qtl_id,
                                        phenotype_info = pheno,
                                        significance_threshold = plot_sig_threshold,
                                        significance_label = significance_label,
                                        title_phenotype_field = title_phenotype_field,
                                        plot_width = plot_width,
                                        plot_height = plot_height,
                                        credible_set = credible_set_snps
                                    )
                                    pname_pdf <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_coloc.pdf"))
                                    pname_png <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_coloc.png"))
                                    pdf_path <- file.path(dir_plots, pname_pdf)
                                    png_path <- file.path(dir_plots, pname_png)
                                    save_res <- save_plot_with_fallback(
                                        plot_obj = p,
                                        pdf_path = pdf_path,
                                        png_path = png_path,
                                        plot_width = plot_width,
                                        plot_height = plot_height
                                    )
                                    plot_path <- save_res$path

                                    if (identical(save_res$format, "png")) {
                                        message(glue("[PLOT] Saved as PNG: {basename(save_res$path)}"))
                                    } else {
                                        message(glue("[PLOT] Saved as PDF: {basename(save_res$path)}"))
                                    }
                                }, error = function(e) {
                                    message(glue("Plot Error: {e$message}"))
                                })

                                tryCatch({
                                    chrom_num <- unique(input_data$CHR.qtl)[1]
                                    susie_best_pp4 <- NA_real_
                                    susie_summary <- NULL
                                    if (!is.null(res$susie_result) &&
                                        !is.null(res$susie_result$summary) &&
                                        nrow(res$susie_result$summary) > 0) {
                                        susie_summary <- res$susie_result$summary
                                        susie_best_pp4 <- max(susie_summary$PP.H4.abf, na.rm = TRUE)
                                    }

                                    rds_data <- list(
                                        merged_data = input_data,
                                        lead_snp = locus$snp,
                                        gene_symbol = gene_sym_for_filename,
                                        gwas_id = gwas_cfg$id,
                                        qtl_id = qtl_id,
                                        chrom = chrom_num,
                                        plink_bfile = plink_ref_hg38,
                                        credible_set = credible_set_snps,
                                        pp4 = pp4,
                                        susie_best_pp4 = susie_best_pp4,
                                        susie_summary = susie_summary
                                    )
                                    rds_fname <- sanitize_filename(glue("{gwas_cfg$id}_{qtl_id}_{gene_sym_for_filename}_coloc.rds"))
                                    rds_path <- file.path(dir_rds, rds_fname)
                                    saveRDS(rds_data, rds_path)
                                    message(glue("[RDS] Saved: {rds_fname}"))
                                }, error = function(e) {
                                    message(glue("[RDS] Save failed: {e$message}"))
                                })
                            }

                            list(
                                record = list(
                                    task_id = task_id,
                                    status = "completed",
                                    gwas_id = gwas_cfg$id,
                                    locus_id = locus_id,
                                    qtl_id = qtl_id,
                                    phenotype = pheno,
                                    identification_method = locus$identification_method,
                                    pp4 = pp4,
                                    n_snps = nsnps,
                                    result_file = results_path,
                                    plot_file = plot_path,
                                    rds_file = rds_path,
                                    error_message = NA_character_
                                ),
                                row = data.table(
                                    GWAS_ID = gwas_cfg$id,
                                    QTL_ID = qtl_id,
                                    Locus = locus$snp,
                                    Phenotype = pheno,
                                    PP4 = pp4,
                                    n_snps = nsnps,
                                    identification_method = locus$identification_method
                                )
                            )
                        }, error = function(e) {
                            list(
                                record = list(
                                    task_id = task_id,
                                    status = "failed",
                                    gwas_id = gwas_cfg$id,
                                    locus_id = locus_id,
                                    qtl_id = qtl_id,
                                    phenotype = pheno,
                                    identification_method = locus$identification_method,
                                    pp4 = NA_real_,
                                    n_snps = NA_real_,
                                    result_file = results_path,
                                    plot_file = NA_character_,
                                    rds_file = NA_character_,
                                    error_message = e$message
                                ),
                                row = NULL
                            )
                        })

                        records[[length(records) + 1]] <- task_result$record
                        if (!is.null(task_result$row)) rows[[length(rows) + 1]] <- task_result$row
                    }

                    list(records = records, rows = rows)
                }

                locus_results <- if (parallel_qtl_tasks && n_cores > 1L && nrow(qtl_meta) >= qtl_parallel_min) {
                    parallel::mclapply(seq_len(nrow(qtl_meta)), process_qtl_wrapper, mc.cores = n_cores)
                } else {
                    lapply(seq_len(nrow(qtl_meta)), process_qtl_wrapper)
                }

                record_list <- unlist(lapply(locus_results, `[[`, "records"), recursive = FALSE)
                row_list <- unlist(lapply(locus_results, `[[`, "rows"), recursive = FALSE)
                final_dt <- if (length(row_list) > 0) rbindlist(row_list, fill = TRUE) else data.table()

                if (nrow(final_dt) > 0) {
                    final_dt <- unique(final_dt)
                    final_dt <- final_dt[order(-PP4)]
                    fwrite(final_dt, results_path)
                    message(glue("-> Saved results for {locus$snp}"))
                }

                if (length(record_list) > 0) {
                    for (rec in record_list) {
                        runtime_update_task(
                            task_id = rec$task_id,
                            status = rec$status,
                            gwas_id = rec$gwas_id,
                            locus_id = rec$locus_id,
                            qtl_id = rec$qtl_id,
                            phenotype = rec$phenotype,
                            identification_method = rec$identification_method,
                            result_file = rec$result_file,
                            plot_file = rec$plot_file,
                            rds_file = rec$rds_file,
                            pp4 = rec$pp4,
                            n_snps = rec$n_snps,
                            error_message = rec$error_message,
                            persist = FALSE
                        )

                        if (write_task_events && rec$status %in% c("completed", "failed", "skipped")) {
                            append_runtime_event(
                                level = if (identical(rec$status, "failed")) "ERROR" else "INFO",
                                stage = paste0("task_", rec$status),
                                message_text = if (identical(rec$status, "failed")) rec$error_message else rec$phenotype,
                                gwas_id = rec$gwas_id,
                                locus_id = rec$locus_id,
                                qtl_id = rec$qtl_id,
                                phenotype = rec$phenotype,
                                task_id = rec$task_id
                            )
                        }
                    }
                    persist_runtime_task_state()
                }

                write_runtime_heartbeat(
                    stage = "locus_complete",
                    message_text = locus_id,
                    counters = list(
                        gwas_index = gwas_idx,
                        total_gwas = total_gwas,
                        locus_index = locus_idx,
                        total_loci = length(loci_list),
                        locus_rows = nrow(final_dt)
                    )
                )
                    }
                }

                append_runtime_event(
                    stage = "gwas_complete",
                    message_text = "Completed GWAS dataset",
                    gwas_id = gwas_cfg$id
                )
            }
        }, error = function(e) {
            warning(glue("[GWAS] {gwas_cfg$id} failed: {e$message}"))
            append_runtime_event(
                level = "ERROR",
                stage = "gwas_failed",
                message_text = e$message,
                gwas_id = gwas_cfg$id
            )
        })
    }

    message("\n==============================================================")
    message("[SUM] Summary all results...")
    message("==============================================================")
    tryCatch({
        merge_all_results(
            output_dir = base_out_dir,
            pp4_threshold = coloc_pp4_thresh,
            merge_susie = TRUE,
            save_summary = TRUE
        )
        append_runtime_event(stage = "summary_complete", message_text = "Merged ABF/SuSiE outputs")
    }, error = function(e) {
        warning(glue("[SUM] Failed to merge results: {e$message}"))
        append_runtime_event(level = "ERROR", stage = "summary_failed", message_text = e$message)
    })

    if (!is.null(cfg_global$sensitivity_analysis) &&
        isTRUE(cfg_global$sensitivity_analysis$enabled)) {
        message("\n==============================================================")
        message("[SENS] Running Prior Sensitivity Analysis...")
        message("==============================================================")
        tryCatch({
            source("src/utils_sensitivity.R")
            abf_files <- list.files(dir_abf, pattern = "_results\\.csv$", full.names = TRUE)
            if (length(abf_files) > 0) {
                example_res <- fread(abf_files[1])
                if (nrow(example_res) > 0 && any(example_res$PP4 >= 0.5, na.rm = TRUE)) {
                    message("[SENS] Sensitivity analysis module loaded")
                    message("[SENS] To run sensitivity analysis, use: run_sensitivity_analysis() on specific loci")
                }
            }
        }, error = function(e) {
            warning(glue("[SENS] Sensitivity analysis failed: {e$message}"))
        })
    }

    message("\n==============================================================")
    message("[REPORT] Generating Interactive HTML Report...")
    message("==============================================================")
    tryCatch({
        source("src/utils_report.R")

        report_file <- file.path(base_out_dir, "coloc_report.html")
        success <- generate_html_report(
            results_dir = base_out_dir,
            output_file = report_file,
            project_name = paste0("EasyColoc Analysis - ", total_gwas, " GWAS datasets")
        )

        if (success) {
            message(glue("[REPORT] Interactive report saved: {report_file}"))
            append_runtime_event(stage = "report_complete", message_text = report_file)
        }
    }, error = function(e) {
        warning(glue("[REPORT] Failed to generate report: {e$message}"))
        append_runtime_event(level = "ERROR", stage = "report_failed", message_text = e$message)
    })

    write_runtime_heartbeat(
        stage = "pipeline_complete",
        message_text = "EasyColoc pipeline finished",
        counters = list(total_gwas = total_gwas)
    )
}
run_pipeline()
message("==============================================================")
message("EasyColoc Analysis Complete!")
message("==============================================================")
