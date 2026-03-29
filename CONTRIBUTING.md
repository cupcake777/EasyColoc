# Contributing

## Development Baseline

Before opening a PR, run:

```bash
Rscript tests/check_parse.R
bash tests/smoke_test_cli.sh
```

If the local environment is fully configured, also run:

```bash
./easycoloc smoke
```

## Contribution Priorities

- Reduce setup friction for new users
- Keep configs portable and explicit
- Prefer machine-readable outputs over log scraping
- Preserve backward compatibility for `Rscript run_coloc.r` when possible

## Review Standard

- Every new entrypoint should have at least a parse or smoke-level check
- Documentation changes should match actual command behavior
- Do not hard-code lab-specific absolute paths into reusable examples
