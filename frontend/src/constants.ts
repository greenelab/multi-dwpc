import type { BiologicalProcess } from "./api/client";

// Worked example: pyrimidine nucleotide catabolic process gene set.
export const EXAMPLE_GENES = [
  "DUT", "ENTPD4", "MBD4", "NEIL1", "NEIL2", "NT5C", "NT5M", "NTHL1", "OGG1",
  "SMUG1", "TDG", "UNG", "DPYD", "DPYS", "TYMP", "UPB1", "UPP1", "UPP2",
].join("\n");

export const EXAMPLE_TARGET: BiologicalProcess = {
  identifier: "GO:0006244",
  name: "pyrimidine nucleotide catabolic process",
};

export const TYPE_COLORS: Record<string, string> = {
  G: "#74c476",
  BP: "#fd8d3c",
  CC: "#fdae6b",
  MF: "#fdd0a2",
  PW: "#9ecae1",
  C: "#c994c7",
  D: "#fb6a4a",
  A: "#bcbddc",
};
export const DEFAULT_NODE_COLOR = "#d9d9d9";

export const NODE_TYPE_NAMES: Record<string, string> = {
  G: "Gene",
  A: "Anatomy",
  BP: "Biological Process",
  CC: "Cellular Component",
  C: "Compound",
  D: "Disease",
  MF: "Molecular Function",
  PC: "Pharmacologic Class",
  PW: "Pathway",
  SE: "Side Effect",
  S: "Symptom",
};
