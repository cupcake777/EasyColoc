# Config Layout

EasyColoc now treats configuration as three layers:

| Layer | Location | Intended use |
| --- | --- | --- |
| Public defaults | `config/*.yaml` | safe, portable repository defaults with no lab-specific absolute paths |
| Local private overrides | `config/local/*.yaml` | user- or lab-specific paths and runtime settings that should not be committed |
| Project template | `templates/project/config/*.yaml` | scaffold copied into a new analysis project via `./easycoloc init` |

## Recommended Practice

1. Leave `config/*.yaml` reusable and publication-safe.
2. Put machine-specific paths in `config/local/`.
3. For real analyses, prefer a standalone project created with `./easycoloc init`.

## Using Local Private Configs

EasyColoc already supports config selection through CLI flags or environment variables:

```bash
./easycoloc run \
  --global config/local/global.yaml \
  --gwas config/local/gwas.yaml \
  --qtl config/local/qtl.yaml
```

or:

```bash
export EASYCOLOC_GLOBAL_CONFIG=config/local/global.yaml
export EASYCOLOC_GWAS_CONFIG=config/local/gwas.yaml
export EASYCOLOC_QTL_CONFIG=config/local/qtl.yaml
./easycoloc run --managed
```

## Notes

- `config/gwas_batch.yaml` remains an example batch-processing file.
- `config/qtl_gtex_generated.yaml` is a generated example produced by bootstrap helpers.
