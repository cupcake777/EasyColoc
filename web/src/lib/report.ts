import type { ReportRow, RuntimeHeartbeat } from "../types";

export function filterResults(
  rows: ReportRow[],
  filters: { query: string; minPp4: number }
): ReportRow[] {
  const needle = filters.query.trim().toLowerCase();

  return rows.filter((row) => {
    const matchesQuery =
      needle.length === 0 ||
      [row.gwas_id, row.qtl_id, row.locus, row.phenotype, row.source_file].some((value) =>
        String(value ?? "").toLowerCase().includes(needle)
      );

    return matchesQuery && row.pp4 >= filters.minPp4;
  });
}

export function topHits(rows: ReportRow[], limit: number): ReportRow[] {
  if (limit <= 0) {
    return [];
  }

  return [...rows].sort((left, right) => right.pp4 - left.pp4).slice(0, limit);
}

export function rowKey(row: ReportRow): string {
  return [row.gwas_id, row.qtl_id, row.locus, row.phenotype].join("::");
}

export function formatPp4(value: number): string {
  return Number.isFinite(value) ? value.toFixed(2) : "NA";
}

export function formatHeartbeatStage(heartbeat: RuntimeHeartbeat | undefined): string | null {
  if (!heartbeat || typeof heartbeat !== "object") {
    return null;
  }

  const stage = heartbeat.stage;
  if (typeof stage === "string" && stage.trim().length > 0) {
    return stage;
  }
  if (Array.isArray(stage)) {
    const parts = stage.filter((value): value is string => typeof value === "string" && value.trim().length > 0);
    return parts.length > 0 ? parts.join(" / ") : null;
  }

  if (heartbeat.missing === true) {
    return "heartbeat missing";
  }
  if (heartbeat.invalid_json === true) {
    return "heartbeat invalid";
  }

  return null;
}

export function formatHeartbeatMessage(heartbeat: RuntimeHeartbeat | undefined): string | null {
  if (!heartbeat || typeof heartbeat !== "object") {
    return null;
  }

  const message = heartbeat.message_text;
  return typeof message === "string" && message.trim().length > 0 ? message : null;
}
