import type { Heatmap as HeatmapData } from "../api/client";

export function Heatmap({ data, metapath }: { data: HeatmapData; metapath: string }) {
  const template = `repeat(${data.genes.length}, 35px)`;
  return (
    <figure className="heatmap">
      <figcaption>Intermediate sharing: {metapath}</figcaption>
      <div className="heatmap-body">
        <div className="heatmap-axis-y">Intermediate</div>
        <div className="heatmap-rows">
          {data.intermediates.map((item) => (
            <div key={item.id} title={item.name}>
              {item.name}
            </div>
          ))}
        </div>
        <div>
          <div
            className="heatmap-grid"
            role="img"
            aria-label={`Genes using each intermediate for ${metapath}`}
            style={{ gridTemplateColumns: template }}
          >
            {data.cells.flatMap((row, rowIndex) =>
              row.map((used, columnIndex) => (
                <div
                  key={`${rowIndex}-${columnIndex}`}
                  className={used ? "heatmap-used" : undefined}
                  title={`${data.genes[columnIndex].name} ${used ? "uses" : "does not use"} ${data.intermediates[rowIndex].name}`}
                />
              )),
            )}
          </div>
          <div className="heatmap-genes" style={{ gridTemplateColumns: template }}>
            {data.genes.map((gene) => (
              <div key={gene.id}>
                <span>{gene.name}</span>
              </div>
            ))}
          </div>
          <div className="heatmap-axis-x">Gene</div>
        </div>
        <div className="heatmap-legend" aria-hidden="true">
          <span>uses intermediate</span>
          <div />
          <span>absent</span>
        </div>
      </div>
    </figure>
  );
}
