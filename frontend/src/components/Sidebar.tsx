import { keepPreviousData, useQuery } from "@tanstack/react-query";
import { useEffect, useState } from "react";
import {
  randomBiologicalProcess,
  randomGenes,
  searchBiologicalProcesses,
  type BiologicalProcess,
} from "../api/client";
import { EXAMPLE_GENES, EXAMPLE_TARGET } from "../constants";
import { Combobox } from "./Combobox";

const SEARCH_DEBOUNCE_MS = 200;

function TargetSelect({
  target,
  onChange,
}: {
  target: BiologicalProcess;
  onChange: (target: BiologicalProcess) => void;
}) {
  const [text, setText] = useState("");
  const [search, setSearch] = useState("");

  useEffect(() => {
    const timer = setTimeout(() => setSearch(text.trim()), SEARCH_DEBOUNCE_MS);
    return () => clearTimeout(timer);
  }, [text]);

  const results = useQuery({
    queryKey: ["biological-processes", search],
    queryFn: ({ signal }) => searchBiologicalProcesses(search, signal),
    placeholderData: keepPreviousData,
  });

  return (
    <Combobox
      label="Target Biological Process"
      value={target.name}
      options={results.data?.results ?? []}
      getKey={(option) => option.identifier}
      getLabel={(option) => option.name}
      onSelect={onChange}
      onSearch={setText}
      loading={results.isFetching}
    />
  );
}

interface SidebarProps {
  target: BiologicalProcess;
  onTargetChange: (target: BiologicalProcess) => void;
  genesRaw: string;
  onGenesChange: (genesRaw: string) => void;
  running: boolean;
  onRun: () => void;
  metapaths: string[] | null;
  metapath: string;
  onMetapathChange: (metapath: string) => void;
}

export function Sidebar({
  target,
  onTargetChange,
  genesRaw,
  onGenesChange,
  running,
  onRun,
  metapaths,
  metapath,
  onMetapathChange,
}: SidebarProps) {
  const [randomError, setRandomError] = useState(false);

  const fill = (load: () => Promise<void>) => () => {
    setRandomError(false);
    load().catch(() => setRandomError(true));
  };

  return (
    <aside className="sidebar">
      <h2>Query</h2>
      <TargetSelect target={target} onChange={onTargetChange} />
      <div className="sidebar-buttons">
        <button
          type="button"
          title="Pyrimidine nucleotide catabolism worked example"
          onClick={() => {
            onTargetChange(EXAMPLE_TARGET);
            onGenesChange(EXAMPLE_GENES);
          }}
        >
          Example query
        </button>
        <button
          type="button"
          onClick={fill(async () => onTargetChange(await randomBiologicalProcess()))}
        >
          Random BP
        </button>
        <button
          type="button"
          onClick={fill(async () => onGenesChange((await randomGenes()).genes.join("\n")))}
        >
          Random genes
        </button>
      </div>
      {randomError && <p className="sidebar-error">Could not load a random selection.</p>}
      <div className="field">
        <label htmlFor="genes">
          Gene symbols (one per line or comma-separated; Entrez IDs also accepted)
        </label>
        <textarea
          id="genes"
          value={genesRaw}
          onChange={(event) => onGenesChange(event.target.value)}
        />
      </div>
      <button
        type="button"
        className="primary"
        disabled={running || !genesRaw.trim()}
        onClick={onRun}
      >
        Run query
      </button>

      {metapaths && (
        <>
          <h2>Drill-down metapath</h2>
          <Combobox
            label="Metapath"
            value={metapath}
            options={metapaths}
            getKey={(option) => option}
            getLabel={(option) => option}
            onSelect={onMetapathChange}
          />
        </>
      )}
    </aside>
  );
}
