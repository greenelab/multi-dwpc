import type { components } from "./schema";

type Schemas = components["schemas"];

export type BiologicalProcess = Schemas["BiologicalProcess"];
export type MetapathRanking = Schemas["MetapathRankingResponse"];
export type MetapathDrilldown = Schemas["MetapathDrilldownResponse"];
export type OverallSubgraphs = Schemas["OverallSubgraphsResponse"];
export type Heatmap = Schemas["Heatmap"];
export type Subgraph = Schemas["Subgraph"];

export interface Query {
  target_id: string;
  genes_raw: string;
}

export class ApiError extends Error {
  constructor(
    message: string,
    readonly status: number,
  ) {
    super(message);
  }
}

function errorMessage(body: unknown): string | null {
  if (typeof body === "string") return body;
  if (Array.isArray(body)) return body.join(" ");
  if (body && typeof body === "object") {
    if ("detail" in body) return String(body.detail);
    return Object.entries(body)
      .map(([field, messages]) => `${field}: ${[messages].flat().join(" ")}`)
      .join(" ");
  }
  return null;
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  let response: Response;
  try {
    response = await fetch(`/api/${path}`, init);
  } catch (cause) {
    if (cause instanceof DOMException && cause.name === "AbortError") throw cause;
    throw new ApiError("Could not reach the API.", 0);
  }
  if (response.ok) return response.json() as Promise<T>;

  const body: unknown = await response.json().catch(() => null);
  throw new ApiError(
    errorMessage(body) ?? `The API returned an error (${response.status}).`,
    response.status,
  );
}

function post<T>(path: string, body: object): Promise<T> {
  return request<T>(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

export function searchBiologicalProcesses(q: string, signal?: AbortSignal) {
  const params = new URLSearchParams({ q });
  return request<Schemas["BiologicalProcessList"]>(
    `reference/biological-processes/?${params}`,
    { signal },
  );
}

export const randomBiologicalProcess = () =>
  request<BiologicalProcess>("reference/biological-processes/random/");

export const randomGenes = () => request<Schemas["RandomGenes"]>("reference/genes/random/");

export const fetchMetapathRanking = (query: Query) =>
  post<MetapathRanking>("queries/metapath-ranking/", query);

export const fetchMetapathDrilldown = (query: Query, metapath: string) =>
  post<MetapathDrilldown>("queries/metapath-drilldown/", { ...query, metapath });

export const fetchOverallSubgraphs = (query: Query) =>
  post<OverallSubgraphs>("queries/overall-subgraphs/", query);
