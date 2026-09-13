#!/usr/bin/env python3
"""Verify the laptop app with the full null bundle: memory, latency, behaviour.

1. In-process: time bundle load and query_metapath_z for the worked example
   (18 genes) and 200 random genes, first and repeat; check determinism and
   the unknown-gene error.
2. The Streamlit app, driven by Playwright: run both queries and an
   unknown-gene query while sampling the server's resident memory.

Writes tables and figures into ``--task-dir``. Requires data/null_bundle.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import threading
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import psutil  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.multi_dwpc_query import query_metapath_z  # noqa: E402
from src.null_bundle import NullBundle  # noqa: E402
from src.query_dwpc import QueryDwpc  # noqa: E402

DATA_DIR = REPO_ROOT / "data"


def _app_constant(name: str):
    """Read a literal constant from app.py without running the Streamlit page."""
    import ast

    for node in ast.parse((REPO_ROOT / "app.py").read_text()).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


EXAMPLE_BP_ID = _app_constant("EXAMPLE_BP_ID")
EXAMPLE_GENES = _app_constant("EXAMPLE_GENES")
PORT = 8599
FULL_PRELOAD_GB = 33.6  # all 52 G->BP matrices as CSR, from the npz headers (decisions.md)
SEED = 20260913


def random_symbols(n: int) -> list[str]:
    genes = pd.read_csv(DATA_DIR / "nodes" / "Gene.tsv", sep="\t")
    return sorted(np.random.default_rng(SEED).choice(genes["name"].astype(str), size=n, replace=False))


def symbols_to_ids(symbols: list[str]) -> list[int]:
    genes = pd.read_csv(DATA_DIR / "nodes" / "Gene.tsv", sep="\t")
    lookup = dict(zip(genes["name"].astype(str), genes["identifier"].astype(int)))
    return [lookup[s] for s in symbols if s in lookup]


def in_process(example: list[str], random200: list[str]) -> pd.DataFrame:
    rows = []
    start = time.perf_counter()
    bundle, query_dwpc = NullBundle(DATA_DIR / "null_bundle", DATA_DIR), QueryDwpc(DATA_DIR)
    rows.append({"step": "load bundle + QueryDwpc", "genes": 0, "seconds": time.perf_counter() - start})
    for label, symbols in [("worked example", example), ("random 200", random200)]:
        ids = symbols_to_ids(symbols)
        frames = []
        for run in ("first", "repeat"):
            start = time.perf_counter()
            frames.append(query_metapath_z(ids, EXAMPLE_BP_ID, bundle=bundle, query_dwpc=query_dwpc))
            rows.append({"step": f"{label} query ({run})", "genes": len(ids), "seconds": time.perf_counter() - start,
                         "metapaths_scored": len(frames[-1])})
        pd.testing.assert_frame_equal(frames[0], frames[1])  # deterministic
    try:
        query_metapath_z([999999999], EXAMPLE_BP_ID, bundle=bundle, query_dwpc=query_dwpc)
        raise AssertionError("unknown genes did not raise")
    except ValueError as exc:
        assert "None of the provided gene IDs" in str(exc)
    return pd.DataFrame(rows)


class MemorySampler(threading.Thread):
    def __init__(self, pid: int):
        super().__init__(daemon=True)
        self.process = psutil.Process(pid)
        self.samples: list[tuple[float, float]] = []
        self.stop = threading.Event()

    def run(self):
        start = time.perf_counter()
        while not self.stop.is_set():
            try:
                procs = [self.process] + self.process.children(recursive=True)
                rss = sum(p.memory_info().rss for p in procs) / 1e9
            except psutil.Error:
                break
            self.samples.append((time.perf_counter() - start, rss))
            time.sleep(0.1)

    def peak(self) -> float:
        return max(r for _, r in self.samples)


def wait_idle(page, timeout_s: float = 600):
    """Wait until Streamlit's running indicator has been gone for 2 s."""
    deadline = time.perf_counter() + timeout_s
    quiet_since = None
    while time.perf_counter() < deadline:
        running = page.locator('[data-testid="stStatusWidget"]').count() > 0
        if running:
            quiet_since = None
        elif quiet_since is None:
            quiet_since = time.perf_counter()
        elif time.perf_counter() - quiet_since > 2:
            return
        time.sleep(0.2)
    raise TimeoutError("app did not become idle")


def drive_app(task_dir: Path, random200: list[str]) -> tuple[pd.DataFrame, MemorySampler]:
    from playwright.sync_api import sync_playwright

    server = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "app.py", "--server.headless", "true",
         "--server.port", str(PORT), "--browser.gatherUsageStats", "false"],
        cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    sampler = MemorySampler(server.pid)
    sampler.start()
    rows = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(viewport={"width": 1500, "height": 1000})
            for _ in range(120):
                try:
                    page.goto(f"http://localhost:{PORT}", timeout=2000)
                    break
                except Exception:
                    time.sleep(0.5)
            page.get_by_role("button", name="Run query").wait_for(timeout=120000)
            wait_idle(page)

            def run_query(label: str, genes_text: str | None, expect: str):
                if genes_text is not None:
                    page.locator("textarea").fill(genes_text)
                    page.keyboard.press("Tab")
                    wait_idle(page)
                start = time.perf_counter()
                page.get_by_role("button", name="Run query").click()
                time.sleep(1.0)  # let the rerun start before checking for idle
                page.get_by_text(expect).first.wait_for(timeout=600000)
                wait_idle(page)
                rows.append({"step": label, "seconds": time.perf_counter() - start,
                             "server_rss_gb_after": sampler.samples[-1][1]})
                page.screenshot(path=str(task_dir / "figures" / f"verify_app_{label.replace(' ', '_')}.png"))

            run_query("worked example", None, "Metapath ranking")
            run_query("random 200", "\n".join(random200), "Metapath ranking")
            run_query("unknown genes", "NOTAGENE1\nNOTAGENE2", "No valid gene symbols")
            browser.close()
    finally:
        sampler.stop.set()
        server.terminate()
        server.wait(timeout=30)
    return pd.DataFrame(rows), sampler


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task-dir", type=Path, default=REPO_ROOT / "docs" / "tasks" / "laptop-null-bundle")
    parser.add_argument("--plot-only", action="store_true", help="Redraw the figure from existing tables")
    args = parser.parse_args()
    if args.plot_only:
        plot(args.task_dir)
        return
    (args.task_dir / "tables").mkdir(parents=True, exist_ok=True)
    (args.task_dir / "figures").mkdir(parents=True, exist_ok=True)
    example = EXAMPLE_GENES.split()
    random200 = random_symbols(200)

    timing = in_process(example, random200)
    timing.to_csv(args.task_dir / "tables" / "verify_query_timing.csv", index=False)
    print(timing.to_string(index=False))

    app_rows, sampler = drive_app(args.task_dir, random200)
    app_rows["server_peak_rss_gb"] = sampler.peak()
    app_rows.to_csv(args.task_dir / "tables" / "verify_app.csv", index=False)
    pd.DataFrame(sampler.samples, columns=["seconds", "server_rss_gb"]).to_csv(
        args.task_dir / "tables" / "verify_app_memory_trace.csv", index=False)
    print(app_rows.to_string(index=False))

    plot(args.task_dir)


def plot(task_dir: Path) -> None:
    """Draw the verify figure from the committed tables."""
    timing = pd.read_csv(task_dir / "tables" / "verify_query_timing.csv")
    app_rows = pd.read_csv(task_dir / "tables" / "verify_app.csv")
    fig, (ax_mem, ax_time) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [1, 1.3]})
    ax_mem.bar(["matrix preload\n(previous app)", "null bundle app\n(measured peak)"],
               [FULL_PRELOAD_GB, app_rows.server_peak_rss_gb.max()], color=["#bbbbbb", "#1f77b4"])
    ax_mem.axhline(24, color="black", ls="--", lw=0.8, label="laptop RAM (24 GB)")
    ax_mem.set_ylabel("resident memory (GB)")
    ax_mem.set_title("Memory, all 52 G->BP metapaths")
    ax_mem.legend(loc="upper right")
    queries = timing[timing.step.str.contains("query")]
    ax_time.barh(queries.step, queries.seconds, color="#1f77b4")
    ax_time.set_xlabel("wall time (s)")
    ax_time.set_title("query_metapath_z, all 52 metapaths")
    fig.tight_layout()
    fig.savefig(task_dir / "figures" / "verify_memory_and_latency.png", dpi=150)


if __name__ == "__main__":
    main()
