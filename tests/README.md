# Example Test Output

This directory contains example output from a test run of EasyColoc.

## Contents

- `abf/` - ABF colocalization results
- `susie/` - SuSiE fine-mapping results
- `plots/` - Visualization plots (LocusZoom-style)

## Usage

These files are provided to show the expected output format. To generate your own results:

1. Configure the files in `config/` with your data paths
2. Run: `Rscript run_coloc.r`
3. Results will be saved to the `output_dir` specified in `config/global.yaml`
