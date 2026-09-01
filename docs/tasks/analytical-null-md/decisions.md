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
