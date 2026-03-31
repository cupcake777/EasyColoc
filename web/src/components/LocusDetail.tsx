import { formatPp4 } from "../lib/report";
import type { ReportRow } from "../types";

type LocusDetailProps = {
  row: ReportRow | null;
};

export function LocusDetail({ row }: LocusDetailProps) {
  if (!row) {
    return (
      <aside className="detail-panel detail-panel--empty">
        <p className="eyebrow">Locus Detail</p>
        <h2>Select a locus</h2>
        <p className="detail-panel__copy">
          Use the table or top-hit strip to inspect one coloc signal in detail.
        </p>
      </aside>
    );
  }

  return (
    <aside className="detail-panel">
      <p className="eyebrow">Locus Detail</p>
      <h2>{row.locus}</h2>
      <p className="detail-panel__copy">
        {row.phenotype} in {row.qtl_id} against {row.gwas_id}
      </p>
      <dl className="detail-grid">
        <div>
          <dt>GWAS</dt>
          <dd>{row.gwas_id}</dd>
        </div>
        <div>
          <dt>QTL</dt>
          <dd>{row.qtl_id}</dd>
        </div>
        <div>
          <dt>Phenotype</dt>
          <dd>{row.phenotype}</dd>
        </div>
        <div>
          <dt>PP4</dt>
          <dd>{formatPp4(row.pp4)}</dd>
        </div>
        <div>
          <dt>n_snps</dt>
          <dd>{row.n_snps}</dd>
        </div>
        <div>
          <dt>Source</dt>
          <dd>{row.source_file || "NA"}</dd>
        </div>
      </dl>
    </aside>
  );
}
