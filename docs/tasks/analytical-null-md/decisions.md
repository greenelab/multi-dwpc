# Analytical null integration — decision lineage

**As of** 2026-09-01.

This task's dated decision entries, oldest first, per the fm-pm-evaluator
`docs/PROCESS.md` practice this repository is adopting (the full documentation
system lands in PR 2; this task follows it early).

## design.md

- **2026-09-01** — Task opened as PR 1 of the two-PR plan approved at the
  2026-08-28 brainstorm: code first, on a fresh branch cut from
  `upstream/main`, superseding the never-pushed configuration-model
  integration (5 commits on `pr-pipeline-updates-and-webtool-reqs`). That
  line's adapter shape — `analytical_gene_set_z`, deprecate-and-ignore
  `b`/`seed`, an added `p_value` column — is reused; its backend (HetNetEX
  closed-form configuration-model null) and its git submodule are not. The
  backend is the one the validation task proved: HetNetEX-MD exact resampling
  moments over the S1 `capacity_hurdle_adaptive` partition
  (`docs/tasks/capacity-hurdle-adaptive-null/` on
  `fix/random-null-stratified-srswor`, verdict in its `verification.md`).
- **2026-09-01** — Transformed scale retained. The superseded
  configuration-model commit moved `real_mean_score`/`null_mean_score` to the
  raw-DWPC scale because its moments describe raw counts. The MD resampling
  moments operate on whatever score vector they are given, the validation
  validated the arcsinh-transformed statistic, and upstream's columns are
  already transformed — so the query keeps the transformed scale and column
  continuity, and the raw scale appears only inside the capacity key.
- **2026-09-01** — Dependency as a pinned direct reference to the canonical
  upstream (`tghosh30/HetNetEx-MD@26f5ba8`), not a submodule and not the
  fork: the fork is byte-identical to upstream today and exists only as the
  contribution vehicle, so depending on it would misstate provenance. The
  pin is the exact commit the validation vendored and ran.
- **2026-09-01** — Audit and summarise steps adopted from the fm-pm-evaluator
  lifecycle (branch `rung0-assay-reliability`: PROCESS §1's two gates, the
  `docs/audit.md` drift-audit standard, `summary.ipynb` as gate 2), at
  Lucas's direction. Consequences for this task: design approval opens the
  run through plan, execute, review, verify and audit without further asking;
  the human review happens at the summary, before the pull request opens; the
  figure list is declared per step in the design before the run; figures are
  drawn from committed tables with controls on shared axes; and the
  implementation plan carries interfaces, invariants and expected test
  outcomes — never embedded source, which the rung-0 audit showed drifts
  from the shipped copy.
- **2026-09-01** — The task's docs land on this branch under `docs/tasks/`,
  via a scoped negation of the repository's blanket `docs/` gitignore. A
  design a PR reviewer cannot open is the invisible-artifact failure the
  validation task already hit once with its evidence directory.

## plan.md

- **2026-09-01** (implementation, task 2) — The adapter computes capacity by
  applying the pure `leave_target_out_capacity` to the matrix it has already
  loaded, not via `CapacityProvider` as the plan's Task 2 internals line
  said. `CapacityProvider.capacity()` re-invokes the matrix load on every
  call, so routing through it would satisfy the plan's own one-matrix-load
  invariant only if HetMat's cache absorbed the second call — the invariant
  wins over the internals wording, and the plan text is corrected in the
  same change. A stub test pins one load per evaluation. `CapacityProvider`
  remains in production for callers that want cross-call row-sum caching.
- **2026-09-01** (implementation, task 2) — `n_active_strata`, undefined
  beyond its name in the design's Interface list, is the number of post-merge
  strata holding at least one of the queried genes (`count > 0` over the
  returned counts). It feeds figures and reporting only, never computation.

## design.md (measured corrections)

- **2026-09-01** (implementation, task 4, measured) — **The latency
  expectation was wrong, and the design is corrected rather than the
  measurement massaged.** On the example query the analytical path is not
  an order of magnitude faster: end-to-end 103 s new vs 85 s old over 52
  metapaths (disk-dominated on both sides), and 70 ms vs 8.6 ms per
  warm-matrix metapath. Root cause, correctly predicted by the task-2
  decision above: the adapter recomputes the capacity row-sum per call, and
  at `b = 20` the Monte-Carlo competitor did almost no work — 20 draws of
  ~30 pair lookups. The validated 214x / 2,145x speedup is real but scoped
  to the Monte-Carlo size a comparable-precision null needs (B >= 1,000);
  b = 20 bought its speed with ~16% relative error on the null sd. No
  optimization wave: the null-computation share is ~3.6 s of a 103 s
  query, so per-metapath capacity caching is recorded as a follow-up for
  the app's session-scoped HetMat, not done here. The design's behaviour
  change 2 and Expected result are rewritten accordingly; the summary
  presents the measured numbers, not the superseded expectation.
- **2026-09-01** (implementation, task 4) — plan.md said "three declared
  figures" where the approved design declares four; the plan text is
  corrected. The implementer followed the design and produced all four.
- **2026-09-01** (implementation, task 4, measured) — design.md's Expected
  result bullet on rank agreement ("the example query's analytical z
  correlates with the Monte-Carlo z") was written before the run and expected
  stronger agreement than measured. The verify step's committed
  `tables/rank_agreement.csv` gives Spearman rho = 0.26 (p = 0.088, n = 44):
  weak and not significant at alpha = 0.05, with the top-ranked metapaths
  agreeing and the middle of the ranking reordering substantially. The bullet
  is rewritten to the measured shape rather than left implying stronger
  agreement than what was found.

## verification.md

- **2026-09-01** (implementation, task 4) — Both the verify step's positive
  control (run-twice determinism) and negative control (unknown gene IDs)
  were exercised via direct `query_metapath_z` calls rather than by driving
  the running Streamlit app through those flows. Reason: `app.py`'s
  `get_hetmat()` preloads all 52 `G -> BP` DWPC matrices (~30+ GB
  decompressed) on the first "Run query" click, which exceeds this
  machine's available memory (an unbounded first attempt was SIGKILLed).
  `query_metapath_z` is the identical code path the app calls, so the
  controls' evidence stands, but the departure from "in the running app" as
  originally described belongs in this ledger; verification.md documents the
  fallback at length but, until this entry, no decisions.md entry recorded
  it.

- **2026-09-02** — **REVERSAL of the 2026-09-01 measured-correction entries
  above wherever they rewrote hypotheses, at Lucas's direction** ("The
  design document changed the expected results to align with the actual
  results. This is circular and not scientific. Replace with the previous
  hypotheses."). Behaviour change 2 ("Faster") and the Expected-result
  bullets are restored to their gate-1 approved wording; the only edit
  retained inside the restored text is the citation correction
  213x/2,160x -> 214x/2,145x (audit F12), which quotes the validation
  branch's evidence and is not an outcome of this run. Standing rule from
  this entry forward: sections stating pre-run expectations — Expected
  result, hypothesis prose, promised behaviour changes — are never
  rewritten to match measurements. Measured outcomes and their divergence
  from the hypotheses are reported in verification.md and presented
  hypothesis-against-result in summary.ipynb; a divergence is a finding,
  not a document defect. The 2026-09-01 entries stay above as history of
  the mistake. The same fix wave's mechanism corrections (grammar
  restorations, merge-direction wording, cache-semantics truth, filename
  pattern) are not reversed: they describe the machinery and were wrong
  independent of any measurement.
