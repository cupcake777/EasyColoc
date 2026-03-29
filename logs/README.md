# Logs

This directory stores local runtime and debug logs that should stay out of the
repository root.

- `managed/`: copied or archived managed-run logs such as `run_easycoloc_full_*.log`
- `pipeline/`: historical pipeline execution logs and troubleshooting reruns
- `plots/`: plot rerun logs and graphics-device byproducts such as `Rplots.pdf`

These files are ignored by Git via the repository `.gitignore`.
