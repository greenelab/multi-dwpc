import { useMutation, useQuery } from "@tanstack/react-query";
import { useId, useState } from "react";
import {
  ApiError,
  fetchMetapathDrilldown,
  fetchMetapathRanking,
  fetchOverallSubgraphs,
  type BiologicalProcess,
  type MetapathRanking,
  type Query,
} from "./api/client";
import { Alert, Spinner } from "./components/Feedback";
import { Sidebar } from "./components/Sidebar";
import { EXAMPLE_GENES, EXAMPLE_TARGET } from "./constants";
import { IntermediateSharingTab } from "./tabs/IntermediateSharingTab";
import { MetapathRankingTab } from "./tabs/MetapathRankingTab";
import { MetapathSubgraphsTab } from "./tabs/MetapathSubgraphsTab";
import { OverallSubgraphsTab } from "./tabs/OverallSubgraphsTab";

const TABS = [
  "Metapath ranking",
  "Intermediate sharing",
  "Overall subgraphs",
  "Metapath subgraphs",
] as const;
type Tab = (typeof TABS)[number];

interface Result {
  query: Query;
  ranking: MetapathRanking;
}

const messageOf = (error: unknown) =>
  error instanceof ApiError ? error.message : "Something went wrong running the query.";

export default function App() {
  const tabsId = useId();
  const [target, setTarget] = useState<BiologicalProcess>(EXAMPLE_TARGET);
  const [genesRaw, setGenesRaw] = useState(EXAMPLE_GENES);
  const [result, setResult] = useState<Result | null>(null);
  const [metapath, setMetapath] = useState("");
  const [tab, setTab] = useState<Tab>(TABS[0]);

  const run = useMutation({
    mutationFn: fetchMetapathRanking,
    onSuccess: (ranking, query) => {
      const metapaths = ranking.metapath_ranking.map((row) => row.metapath);
      setResult({ query, ranking });
      setMetapath((current) => (metapaths.includes(current) ? current : metapaths[0]));
    },
    onError: () => setResult(null),
  });

  const query = result?.query;
  const drilldown = useQuery({
    queryKey: ["metapath-drilldown", query, metapath],
    queryFn: () => fetchMetapathDrilldown(query!, metapath),
    enabled: Boolean(query && metapath),
  });
  const overall = useQuery({
    queryKey: ["overall-subgraphs", query],
    queryFn: () => fetchOverallSubgraphs(query!),
    enabled: Boolean(query) && tab === "Overall subgraphs",
  });

  return (
    <div className="app">
      <Sidebar
        target={target}
        onTargetChange={setTarget}
        genesRaw={genesRaw}
        onGenesChange={setGenesRaw}
        running={run.isPending}
        onRun={() => run.mutate({ target_id: target.identifier, genes_raw: genesRaw })}
        metapaths={result?.ranking.metapath_ranking.map((row) => row.metapath) ?? null}
        metapath={metapath}
        onMetapathChange={setMetapath}
      />
      <main>
        <h1>Multi-DWPC query tool</h1>
        <p className="caption">Gene list → target → metapath z, intermediate sharing, surviving subpaths</p>

        {run.isPending && <Spinner>Computing metapath z (analytical null)</Spinner>}
        {run.isError && <Alert kind="error">{messageOf(run.error)}</Alert>}
        {!result && !run.isPending && !run.isError && (
          <Alert kind="info">Enter a target and gene list in the sidebar, then run the query.</Alert>
        )}

        {result && (
          <>
            <p className="caption">Last query runtime: {result.ranking.elapsed_s.toFixed(2)} seconds</p>
            <div className="tabs" role="tablist">
              {TABS.map((name) => (
                <button
                  key={name}
                  type="button"
                  role="tab"
                  id={`${tabsId}-${name}`}
                  aria-selected={tab === name}
                  aria-controls={`${tabsId}-${name}-panel`}
                  onClick={() => setTab(name)}
                >
                  {name}
                </button>
              ))}
            </div>
            {TABS.map((name) => (
              <section
                key={name}
                role="tabpanel"
                id={`${tabsId}-${name}-panel`}
                aria-labelledby={`${tabsId}-${name}`}
                hidden={tab !== name}
              >
                {name === "Metapath ranking" && (
                  <MetapathRankingTab rows={result.ranking.metapath_ranking} />
                )}
                {name === "Overall subgraphs" && (
                  <>
                    {overall.isPending && overall.fetchStatus !== "idle" && (
                      <Spinner>Pooling paths across significant metapaths</Spinner>
                    )}
                    {overall.isError && <Alert kind="error">{messageOf(overall.error)}</Alert>}
                    {overall.data && <OverallSubgraphsTab data={overall.data} />}
                  </>
                )}
                {(name === "Intermediate sharing" || name === "Metapath subgraphs") && (
                  <>
                    {drilldown.isPending && <Spinner>Enumerating paths for {metapath}</Spinner>}
                    {drilldown.isError && <Alert kind="error">{messageOf(drilldown.error)}</Alert>}
                    {drilldown.data &&
                      (name === "Intermediate sharing" ? (
                        <IntermediateSharingTab data={drilldown.data} />
                      ) : (
                        <MetapathSubgraphsTab data={drilldown.data} />
                      ))}
                  </>
                )}
              </section>
            ))}
          </>
        )}
      </main>
    </div>
  );
}
