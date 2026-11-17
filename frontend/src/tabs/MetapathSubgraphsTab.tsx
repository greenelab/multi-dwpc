import type { MetapathDrilldown } from "../api/client";
import { DataTable, type Column } from "../components/DataTable";
import { Alert, Expander } from "../components/Feedback";
import { SubgraphChart } from "../components/SubgraphChart";

function pathColumns(paths: MetapathDrilldown["paths"]): Column[] {
  const hops = new Set<string>();
  for (const path of paths) {
    for (const key of Object.keys(path)) {
      const hop = key.match(/^(hop_\d+)_name$/)?.[1];
      if (hop) hops.add(hop);
    }
  }
  return [
    ...["gene", "metapath", "path_rank", "path_score"].map((key) => ({ key, label: key })),
    ...[...hops].sort().map((hop) => ({ key: `${hop}_name`, label: hop })),
  ];
}

export function MetapathSubgraphsTab({ data }: { data: MetapathDrilldown }) {
  if (!data.paths.length) {
    return (
      <>
        <Alert kind="warning">
          No paths survived the path-z filter (≥ {data.path_z_min.toFixed(2)}).
        </Alert>
        <Expander title="Diagnostics">
          <pre className="json">{JSON.stringify(data.diagnostics, null, 2)}</pre>
        </Expander>
      </>
    );
  }
  return (
    <>
      <h3>Top paths by score: {data.metapath}</h3>
      <p className="caption">
        Node size = summed path score through that intermediate. Shared-intermediate information
        for this metapath is shown as a heatmap in the <em>Intermediate sharing</em> tab.
      </p>
      {data.subgraph ? (
        <SubgraphChart data={data.subgraph} />
      ) : (
        <Alert kind="info">No paths have path_score populated.</Alert>
      )}
      <Expander title="Table of surviving subpaths">
        <DataTable
          name={`surviving-subpaths-${data.metapath}`}
          columns={pathColumns(data.paths)}
          rows={data.paths}
        />
      </Expander>
    </>
  );
}
