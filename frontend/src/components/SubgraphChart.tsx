import { useId } from "react";
import type { Subgraph } from "../api/client";
import { DEFAULT_NODE_COLOR, NODE_TYPE_NAMES, TYPE_COLORS } from "../constants";

const UNIT_X = 240;
const UNIT_Y = 55;
const BOX_WIDTH = 0.9 * UNIT_X;
const BOX_HEIGHT = 0.44 * UNIT_Y;
const TITLE_HEIGHT = 44;
const AXIS_HEIGHT = 40;
const LEGEND_WIDTH = 180;

const colorFor = (type: string) => TYPE_COLORS[type] ?? DEFAULT_NODE_COLOR;

/** Layered graph: hop positions as columns, with nodes stacked and centred in each. */
export function SubgraphChart({ data }: { data: Subgraph }) {
  const markerId = `arrow-${useId().replaceAll(":", "")}`;
  const columnCount = data.columns.length;

  const counts = new Array<number>(columnCount).fill(0);
  for (const node of data.nodes) counts[node.column] += 1;
  const tallest = Math.max(1, ...counts);

  const plotWidth = columnCount * UNIT_X;
  const plotHeight = (tallest + 1.4) * UNIT_Y;
  const width = plotWidth + LEGEND_WIDTH;
  const height = TITLE_HEIGHT + plotHeight + AXIS_HEIGHT;
  const columnX = (column: number) => (column + 0.5) * UNIT_X;
  const middleY = TITLE_HEIGHT + plotHeight / 2;

  const seen = new Array<number>(columnCount).fill(0);
  const positions = new Map(
    data.nodes.map((node) => {
      const row = seen[node.column]++;
      const y = middleY + (row - (counts[node.column] - 1) / 2) * UNIT_Y;
      return [node.id, { x: columnX(node.column), y, halfWidth: (BOX_WIDTH * node.scale) / 2 }];
    }),
  );
  const types = [...new Set(data.nodes.map((node) => node.type))].sort();

  return (
    <figure className="subgraph">
      <svg viewBox={`0 0 ${width} ${height}`} style={{ maxWidth: width }} role="img" aria-label={data.title}>
        <defs>
          <marker id={markerId} markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto">
            <path d="M0 0 7 3.5 0 7" fill="none" stroke="gray" strokeWidth="1" />
          </marker>
        </defs>
        <text className="subgraph-title" x={plotWidth / 2} y={26} textAnchor="middle">
          {data.title}
        </text>

        {data.edges.map((edge) => {
          const source = positions.get(edge.source);
          const target = positions.get(edge.target);
          if (!source || !target) return null;
          return (
            <line
              key={`${edge.source}>${edge.target}`}
              x1={source.x + source.halfWidth}
              y1={source.y}
              x2={target.x - target.halfWidth}
              y2={target.y}
              stroke="gray"
              strokeOpacity="0.5"
              markerEnd={`url(#${markerId})`}
            />
          );
        })}

        {data.nodes.map((node) => {
          const position = positions.get(node.id)!;
          const boxWidth = BOX_WIDTH * node.scale;
          const boxHeight = BOX_HEIGHT * node.scale;
          return (
            <g key={node.id}>
              <title>{node.label}</title>
              <rect
                x={position.x - boxWidth / 2}
                y={position.y - boxHeight / 2}
                width={boxWidth}
                height={boxHeight}
                rx="3"
                fill={colorFor(node.type)}
                stroke="black"
                strokeWidth="0.75"
              />
              <text x={position.x} y={position.y} textAnchor="middle" dominantBaseline="central">
                {node.label}
              </text>
            </g>
          );
        })}

        <line
          x1={0}
          x2={plotWidth}
          y1={height - AXIS_HEIGHT}
          y2={height - AXIS_HEIGHT}
          stroke="black"
          strokeWidth="0.75"
        />
        {data.columns.map((label, column) => (
          <text
            key={column}
            className="subgraph-axis"
            x={columnX(column)}
            y={height - AXIS_HEIGHT + 22}
            textAnchor="middle"
          >
            {label}
          </text>
        ))}

        <g transform={`translate(${plotWidth + 16} ${middleY - (types.length * 20 + 24) / 2})`}>
          <text className="subgraph-legend-title" y={12}>
            Node type
          </text>
          {types.map((type, index) => (
            <g key={type} transform={`translate(0 ${24 + index * 20})`}>
              <rect width="20" height="11" fill={colorFor(type)} stroke="black" strokeWidth="0.75" />
              <text x={28} y={6} dominantBaseline="central">
                {NODE_TYPE_NAMES[type] ?? type}
              </text>
            </g>
          ))}
        </g>
      </svg>
    </figure>
  );
}
