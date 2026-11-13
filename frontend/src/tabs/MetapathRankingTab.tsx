import type { MetapathRanking } from "../api/client";
import { DataTable, type Column } from "../components/DataTable";

const COLUMNS: Column[] = [
  "metapath",
  "real_mean_score",
  "null_mean_score",
  "null_std_score",
  "diff",
  "effect_size_z",
  "p_value",
].map((key) => ({ key, label: key }));

export function MetapathRankingTab({ rows }: { rows: MetapathRanking["metapath_ranking"] }) {
  return <DataTable name="metapath-ranking" columns={COLUMNS} rows={rows} />;
}
