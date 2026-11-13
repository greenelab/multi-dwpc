import type { OverallSubgraphs } from "../api/client";
import { Alert } from "../components/Feedback";
import { SubgraphChart } from "../components/SubgraphChart";

export function OverallSubgraphsTab({ data }: { data: OverallSubgraphs }) {
  const { pool_z_min: zMin, pool_z_used: zUsed, pooled_metapaths: metapaths } = data;
  if (!data.top_shared_subgraph && !data.top_paths_subgraph) {
    return (
      <Alert kind="warning">
        No paths across metapaths with z ≥ {zMin.toFixed(2)} (nor z &gt; 0).
      </Alert>
    );
  }
  const fallback = zUsed < zMin ? ` (fallback: nothing cleared z≥${zMin})` : "";
  return (
    <>
      <p className="caption">
        Pooled over {metapaths.length} metapaths at z ≥ {zUsed.toFixed(2)}
        {fallback}: {metapaths.join(", ")}
      </p>
      <h3>Top shared intermediates</h3>
      <p className="caption">Node size = number of genes using that intermediate.</p>
      {data.top_shared_subgraph ? (
        <SubgraphChart data={data.top_shared_subgraph} />
      ) : (
        <Alert kind="info">No paths survived the top-shared filter.</Alert>
      )}
      <h3>Top paths by score</h3>
      <p className="caption">Node size = summed path score through that intermediate.</p>
      {data.top_paths_subgraph ? (
        <SubgraphChart data={data.top_paths_subgraph} />
      ) : (
        <Alert kind="info">No paths have path_score populated.</Alert>
      )}
    </>
  );
}
