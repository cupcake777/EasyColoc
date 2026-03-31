import { describe, expect, it } from "vitest";
import { filterResults, topHits } from "./report";

const rows = [
  {
    gwas_id: "GWAS_A",
    qtl_id: "postnatal",
    locus: "rs100",
    phenotype: "GENE1",
    pp4: 0.92,
    n_snps: 41
  },
  {
    gwas_id: "GWAS_B",
    qtl_id: "prenatal",
    locus: "rs101",
    phenotype: "GENE2",
    pp4: 0.41,
    n_snps: 25
  }
];

describe("report helpers", () => {
  it("filters by free-text query and minimum pp4", () => {
    const filtered = filterResults(rows, { query: "gene1", minPp4: 0.8 });

    expect(filtered).toHaveLength(1);
    expect(filtered[0].locus).toBe("rs100");
  });

  it("returns top hits in descending PP4 order", () => {
    const hits = topHits(rows, 1);

    expect(hits).toHaveLength(1);
    expect(hits[0].pp4).toBe(0.92);
  });
});
