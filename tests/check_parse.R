#!/usr/bin/env Rscript

files_to_check <- c(
  "run_coloc.r",
  "tools/rerun_plots.R",
  "src/utils_config.R",
  "src/utils_download.R",
  "src/utils_reference.R",
  "src/utils_bootstrap.R",
  "src/utils_output.R",
  "src/utils_runtime.R",
  "src/utils_helpers.R",
  "src/utils_analysis.R",
  "src/utils_report.R",
  "src/utils_sensitivity.R",
  "src/utils_plot.R",
  "tools/doctor_easycoloc.R",
  "tools/bootstrap_references.R",
  "tools/list_reference_requirements.R",
  "tools/monitor_easycoloc.R",
  "tools/check_run_completion.R",
  "tools/build_output_manifest.R",
  "tools/summarize_run_status.R",
  "tests/smoke_test_plotting.R",
  "tests/smoke_test_bootstrap_refs.R",
  "tests/smoke_test_bootstrap_demo.R",
  "tests/smoke_test_gtex_bootstrap.R",
  "tests/smoke_test_1kg_setup.R",
  "tests/smoke_test_ld_cache.R",
  "tests/smoke_test_doctor.R",
  "tests/smoke_test_refs.R",
  "tests/smoke_test_runtime_monitor.R",
  "tests/smoke_test_output_dir_tools.R",
  "examples/minimal/synthetic_plot_demo.R"
)

for (file in files_to_check) {
  parse(file = file)
}

cat("[PARSE] checked", length(files_to_check), "files\n")
