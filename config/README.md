# Config Layout

EasyColoc keeps repository configuration small and portable:

| Layer | Location | Intended use |
| --- | --- | --- |
| Public defaults | `config/global.yml`, `config/gwas.yml`, `config/qtl.yml` | safe, portable repository defaults with no lab-specific absolute paths |
| Local private overrides | `config/local/*.yml` | user- or lab-specific paths and runtime settings that should not be committed |

## Recommended Practice

1. Leave the three public config files reusable and publication-safe.
2. Put machine-specific paths in `config/local/`.
3. For real analyses, create a minimal standalone project with `./easycoloc init`.

## Using Local Private Configs

EasyColoc already supports config selection through CLI flags or environment variables:

```bash
./easycoloc run \
  --global config/local/global.yml \
  --gwas config/local/gwas.yml \
  --qtl config/local/qtl.yml
```

or:

```bash
export EASYCOLOC_GLOBAL_CONFIG=config/local/global.yml
export EASYCOLOC_GWAS_CONFIG=config/local/gwas.yml
export EASYCOLOC_QTL_CONFIG=config/local/qtl.yml
./easycoloc run --managed
```

Generated QTL config files from GTEx metadata should be written to a working
project or local output directory, not committed as repository defaults.
