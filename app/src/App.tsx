import "@fontsource/lato/400.css";
import "@fontsource/lato/700.css";
import "@fontsource-variable/roboto-mono/wght.css";
import "./styles.css";
import { useState } from "react";

const example = `Search 1\nSearch 2\nSearch 3\nSearch 4\nSearch 5`;

export default function App() {
  const [search, setSearch] = useState("");

  return (
    <>
      <header className="flex items-center justify-center bg-gray-800 p-8 text-white">
        <h1 className="text-white">
          <svg viewBox="-50 -50 100 100" className="h-20">
            <g
              fill="none"
              stroke="white"
              strokeWidth="4"
              strokeDasharray="0.0001 8"
              strokeDashoffset="1"
              strokeLinecap="round"
            >
              <path d="M -15 -30 L 0 -30 A 30 30 0 0 1 30 0" />
              <path d="M -40 0 L 30 0" />
              <path d="M -15 30 L 0 30 A 30 30 0 0 0 30 0" />
            </g>
            <circle cx="30" cy="0" r="10" fill="#02b3e4"></circle>
            <circle cx="-15" cy="-30" r="10" fill="#fa700f"></circle>
            <circle cx="-40" cy="0" r="10" fill="#e91e63"></circle>
            <circle cx="-15" cy="30" r="10" fill="#c341d8"></circle>
          </svg>
          Multi-DWPC
        </h1>
      </header>

      <main>
        <section>
          <h2>Source Nodes</h2>

          <div className="grid grid-cols-2 gap-8 max-md:grid-cols-1">
            <textarea
              placeholder={`Search for nodes\nOne per line`}
              rows={5}
              value={search}
              onChange={(event) => setSearch(event.target.value)}
            />
            <div className="border border-gray-200 p-4 leading-loose trim">
              {search
                .trim()
                .split("\n")
                .filter(Boolean)
                .map((line, index) => (
                  <div key={index}>{line} result</div>
                ))}
              {!search.trim() && (
                <div className="text-gray-500">
                  Results
                  <br />
                  will
                  <br />
                  show
                  <br />
                  here
                </div>
              )}
            </div>
          </div>

          <button className="self-center" onClick={() => setSearch(example)}>
            Example
          </button>
        </section>

        <section>
          <h2>Target Node</h2>

          <input placeholder="Search for a node" />
        </section>

        <section>
          <h2>Results</h2>

          <h3>Ranked Multi-DWPC</h3>

          <table>
            <thead>
              <tr>
                <th>Metapath</th>
                <th>Real mean</th>
                <th>Null mean</th>
                <th>Diff</th>
                <th>Z-score</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>Lorem ipsum</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>1</td>
              </tr>
              <tr>
                <td>Lorem ipsum</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>1</td>
              </tr>
              <tr>
                <td>Lorem ipsum</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>0.12345</td>
                <td>1</td>
              </tr>
            </tbody>
          </table>

          <h3>Intermediate Nodes</h3>

          <svg className="size-100 bg-gray-500"></svg>

          <h3>Subgraph</h3>

          <svg className="size-100 bg-gray-500"></svg>
        </section>
      </main>

      <footer className="bg-gray-800 p-8 text-center text-white">
        Project of the Greene Lab
      </footer>
    </>
  );
}
