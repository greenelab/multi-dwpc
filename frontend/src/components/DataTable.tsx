import { useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";

export type Row = Record<string, unknown>;

export interface Column {
  key: string;
  label: string;
}

interface DataTableProps {
  name: string;
  columns: Column[];
  rows: Row[];
}

type Sort = { key: string; direction: 1 | -1 } | null;

const ROW_HEIGHT = 35;
const HEADER_HEIGHT = 36;
const MAX_HEIGHT = 400;
const OVERSCAN = 8;
const INDEX_KEY = "__index";

export function formatCell(value: unknown): string {
  if (value === null || value === undefined) return "";
  if (typeof value !== "number") return String(value);
  if (Number.isInteger(value)) return String(value);
  if (Math.abs(value) < 1e-4) return value.toExponential(2);
  return String(Number(value.toFixed(4)));
}

function compare(a: unknown, b: unknown): number {
  if (typeof a === "number" && typeof b === "number") return a - b;
  return String(a).localeCompare(String(b));
}

function downloadCsv(name: string, columns: Column[], rows: Row[]) {
  const escape = (value: unknown) => {
    const text = value === null || value === undefined ? "" : String(value);
    return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  };
  const lines = [
    columns.map((column) => escape(column.label)).join(","),
    ...rows.map((row) => columns.map((column) => escape(row[column.key])).join(",")),
  ];
  const url = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/csv" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = `${name}.csv`;
  link.click();
  URL.revokeObjectURL(url);
}

function columnWidth(label: string, values: string[]): number {
  const chars = Math.max(label.length + 2, ...values.map((value) => value.length));
  return Math.min(Math.max(chars * 7.5 + 24, 60), 640);
}

function ToolbarIcon({ path }: { path: string }) {
  return (
    <svg viewBox="0 0 24 24" aria-hidden="true">
      <path d={path} />
    </svg>
  );
}

const ICONS = {
  columns: "M3 12s3.5-6 9-6 9 6 9 6-3.5 6-9 6-9-6-9-6Zm9 3a3 3 0 1 0 0-6 3 3 0 0 0 0 6Z",
  download: "M12 3v11m0 0 4-4m-4 4-4-4M5 19h14",
  search: "m20 20-4.5-4.5M17 10.5a6.5 6.5 0 1 1-13 0 6.5 6.5 0 0 1 13 0Z",
  fullscreen: "M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5",
  close: "M6 6l12 12M18 6 6 18",
};

/** Streamlit-style data grid: sorting, search, column toggles, CSV and fullscreen. */
export function DataTable({ name, columns, rows }: DataTableProps) {
  const scrollRef = useRef<HTMLDivElement>(null);
  const [viewHeight, setViewHeight] = useState(MAX_HEIGHT);
  const [scrollTop, setScrollTop] = useState(0);
  const [sort, setSort] = useState<Sort>(null);
  const [hidden, setHidden] = useState<Set<string>>(new Set());
  const [search, setSearch] = useState("");
  const [searchOpen, setSearchOpen] = useState(false);
  const [columnsOpen, setColumnsOpen] = useState(false);
  const [fullscreen, setFullscreen] = useState(false);

  const shownColumns = useMemo(
    () => columns.filter((column) => !hidden.has(column.key)),
    [columns, hidden],
  );

  const displayRows = useMemo(() => {
    const needle = search.trim().toLowerCase();
    let result = rows.map((row, index) => ({ ...row, [INDEX_KEY]: index }) as Row);
    if (needle) {
      result = result.filter((row) =>
        shownColumns.some((column) => formatCell(row[column.key]).toLowerCase().includes(needle)),
      );
    }
    if (sort) {
      const isNull = (row: Row) => row[sort.key] === null || row[sort.key] === undefined;
      result = [...result].sort((a, b) => {
        if (isNull(a) || isNull(b)) return Number(isNull(a)) - Number(isNull(b));
        return compare(a[sort.key], b[sort.key]) * sort.direction;
      });
    }
    return result;
  }, [rows, shownColumns, search, sort]);

  const widths = useMemo(() => {
    const sample = rows.slice(0, 300);
    return shownColumns.map((column) =>
      columnWidth(column.label, sample.map((row) => formatCell(row[column.key]))),
    );
  }, [rows, shownColumns]);

  const indexWidth = String(rows.length).length * 9 + 26;
  const allWidths = [indexWidth, ...widths];
  const template = allWidths.map((width) => `minmax(${width}px, ${width}fr)`).join(" ");
  const minWidth = allWidths.reduce((total, width) => total + width, 0);

  useLayoutEffect(() => {
    const element = scrollRef.current;
    if (!element) return;
    const observer = new ResizeObserver(() => setViewHeight(element.clientHeight));
    observer.observe(element);
    return () => observer.disconnect();
  }, []);

  useEffect(() => {
    if (!fullscreen) return;
    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === "Escape") setFullscreen(false);
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [fullscreen]);

  const toggleSort = (key: string) =>
    setSort((current) => {
      if (current?.key !== key) return { key, direction: 1 };
      return current.direction === 1 ? { key, direction: -1 } : null;
    });

  const toggleColumn = (key: string) =>
    setHidden((current) => {
      const next = new Set(current);
      if (!next.delete(key)) next.add(key);
      return next;
    });

  const first = Math.max(Math.floor(scrollTop / ROW_HEIGHT) - OVERSCAN, 0);
  const last = Math.min(
    Math.ceil((scrollTop + viewHeight) / ROW_HEIGHT) + OVERSCAN,
    displayRows.length,
  );
  const contentHeight = HEADER_HEIGHT + displayRows.length * ROW_HEIGHT + 2;

  const pinned = searchOpen || columnsOpen || fullscreen;
  return (
    <div className={`table${fullscreen ? " table-fullscreen" : ""}${pinned ? " table-pinned" : ""}`}>
      <div className="table-toolbar">
        {searchOpen && (
          <input
            autoFocus
            aria-label={`Search ${name}`}
            placeholder="Search"
            value={search}
            onChange={(event) => setSearch(event.target.value)}
          />
        )}
        <button
          type="button"
          title="Show/hide columns"
          aria-label="Show/hide columns"
          onClick={() => setColumnsOpen((open) => !open)}
        >
          <ToolbarIcon path={ICONS.columns} />
        </button>
        <button
          type="button"
          title="Download as CSV"
          aria-label="Download as CSV"
          onClick={() => downloadCsv(name, shownColumns, displayRows)}
        >
          <ToolbarIcon path={ICONS.download} />
        </button>
        <button
          type="button"
          title="Search"
          aria-label="Search"
          onClick={() => {
            setSearchOpen((open) => !open);
            setSearch("");
          }}
        >
          <ToolbarIcon path={ICONS.search} />
        </button>
        <button
          type="button"
          title={fullscreen ? "Close" : "Fullscreen"}
          aria-label={fullscreen ? "Close" : "Fullscreen"}
          onClick={() => setFullscreen((open) => !open)}
        >
          <ToolbarIcon path={fullscreen ? ICONS.close : ICONS.fullscreen} />
        </button>
        {columnsOpen && (
          <div className="table-columns">
            {columns.map((column) => (
              <label key={column.key}>
                <input
                  type="checkbox"
                  checked={!hidden.has(column.key)}
                  onChange={() => toggleColumn(column.key)}
                />
                {column.label}
              </label>
            ))}
          </div>
        )}
      </div>
      <div
        ref={scrollRef}
        className="table-scroll"
        style={fullscreen ? undefined : { height: Math.min(contentHeight, MAX_HEIGHT) }}
        onScroll={(event) => setScrollTop(event.currentTarget.scrollTop)}
      >
        <div
          role="table"
          aria-label={name}
          aria-rowcount={displayRows.length + 1}
          className="table-grid"
          style={{ minWidth, height: contentHeight - 2 }}
        >
          <div role="row" className="table-row table-header" style={{ gridTemplateColumns: template }}>
            <div role="columnheader" />
            {shownColumns.map((column) => (
              <button
                type="button"
                role="columnheader"
                key={column.key}
                aria-sort={
                  sort?.key !== column.key
                    ? "none"
                    : sort.direction === 1
                      ? "ascending"
                      : "descending"
                }
                onClick={() => toggleSort(column.key)}
              >
                {column.label}
                {sort?.key === column.key && (sort.direction === 1 ? " ↑" : " ↓")}
              </button>
            ))}
          </div>
          {displayRows.slice(first, last).map((row, offset) => {
            const position = first + offset;
            return (
              <div
                role="row"
                className="table-row"
                key={row[INDEX_KEY] as number}
                style={{
                  gridTemplateColumns: template,
                  transform: `translateY(${HEADER_HEIGHT + position * ROW_HEIGHT}px)`,
                }}
              >
                <div role="cell" className="table-index">
                  {row[INDEX_KEY] as number}
                </div>
                {shownColumns.map((column) => {
                  const value = row[column.key];
                  return (
                    <div
                      role="cell"
                      key={column.key}
                      className={typeof value === "number" ? "table-number" : undefined}
                      title={formatCell(value)}
                    >
                      {formatCell(value)}
                    </div>
                  );
                })}
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
}
