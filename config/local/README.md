# Local Private Configs

Use this directory for machine-specific or lab-specific configs that should not
be committed.

Suggested files:

- `config/local/global.yml`
- `config/local/gwas.yml`
- `config/local/qtl.yml`

Invoke them explicitly:

```bash
./easycoloc run \
  --global config/local/global.yml \
  --gwas config/local/gwas.yml \
  --qtl config/local/qtl.yml
```

Or export:

```bash
export EASYCOLOC_GLOBAL_CONFIG=config/local/global.yml
export EASYCOLOC_GWAS_CONFIG=config/local/gwas.yml
export EASYCOLOC_QTL_CONFIG=config/local/qtl.yml
```
