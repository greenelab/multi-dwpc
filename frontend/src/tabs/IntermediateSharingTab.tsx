import type { MetapathDrilldown } from "../api/client";
import { DataTable, type Column } from "../components/DataTable";
import { Alert, Expander } from "../components/Feedback";
import { Heatmap } from "../components/Heatmap";

const COLUMNS: Column[] = ["intermediate", "n_genes_sharing", "genes"].map((key) => ({
  key,
  label: key,
}));

export function IntermediateSharingTab({ data }: { data: MetapathDrilldown }) {
  if (!data.paths.length) {
    return <Alert kind="warning">No intermediates survived the path-z filter.</Alert>;
  }
  if (!data.intermediate_sharing.length) {
    return <Alert kind="info">No intermediate hops for this metapath.</Alert>;
  }
  return (
    <>
      {data.heatmap && <Heatmap data={data.heatmap} metapath={data.metapath} />}
      <Expander title="Intermediate sharing table">
        <DataTable
          name={`intermediate-sharing-${data.metapath}`}
          columns={COLUMNS}
          rows={data.intermediate_sharing}
        />
      </Expander>
    </>
  );
}
