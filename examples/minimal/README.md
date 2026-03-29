# Minimal Example

This directory provides a self-contained plotting demo that does not require
running the full colocalization pipeline.

## What It Shows

- Synthetic locus-level QTL summary data
- A lead SNP plus a small 95% credible set
- The current compact publication-style locus plot
- Automatic subtitle annotation when the phenotype transcript lies outside the
  lead-SNP-centered plotting window

## Run

From the repository root:

```bash
Rscript examples/minimal/synthetic_plot_demo.R
```

Outputs will be written to:

```text
examples/minimal/output/
```

## Why This Exists

This example is meant to be the smallest reproducible entrypoint for:

- checking plotting dependencies
- verifying style changes
- demonstrating the expected figure structure to new users
