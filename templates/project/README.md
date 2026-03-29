# EasyColoc Project Template

This scaffold is meant to be copied outside the EasyColoc source tree and run
via:

```bash
./easycoloc doctor --global /path/to/project/config/global.yaml --gwas /path/to/project/config/gwas.yaml --qtl /path/to/project/config/qtl.yaml
./easycoloc run --managed --global /path/to/project/config/global.yaml --gwas /path/to/project/config/gwas.yaml --qtl /path/to/project/config/qtl.yaml
```

What to edit first:

1. `config/global.yaml`
2. `config/gwas.yaml`
3. `config/qtl.yaml`
4. `config/reference_sources.template.yaml`
5. `data/qtl/QTL_summary.template.csv`

Important:

- This template sets `project_root: ".."` in `config/global.yaml`, so relative
  paths resolve from the project root rather than the `config/` directory.
- To materialize shared references into `refs/`, fill
  `config/reference_sources.template.yaml` and run:
  `./easycoloc bootstrap-refs config/reference_sources.template.yaml --rewrite-config`
- Replace placeholder reference paths before running real analyses.
- Keep real outputs under `results/` and temporary files under `temp/`.
