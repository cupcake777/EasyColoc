# Local Private Configs

Use this directory for machine-specific or lab-specific configs that should not
be committed.

Suggested files:

- `config/local/global.yaml`
- `config/local/gwas.yaml`
- `config/local/qtl.yaml`

Invoke them explicitly:

```bash
./easycoloc run \
  --global config/local/global.yaml \
  --gwas config/local/gwas.yaml \
  --qtl config/local/qtl.yaml
```

Or export:

```bash
export EASYCOLOC_GLOBAL_CONFIG=config/local/global.yaml
export EASYCOLOC_GWAS_CONFIG=config/local/gwas.yaml
export EASYCOLOC_QTL_CONFIG=config/local/qtl.yaml
```
