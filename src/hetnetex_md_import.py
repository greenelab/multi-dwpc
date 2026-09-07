"""Single, monitored import point for HetNetEX-MD (pinned via pyproject.toml
to `hetnetex-md @ git+https://github.com/tghosh30/HetNetEx-MD.git@
26f5ba85bd9b9886eafebcbc544879e531f5262e`).

Per docs/superpowers/specs/2026-08-07-hetnetex-md-validation-design.md's
independent-evaluation verdict, ONLY `exact_resampling_moments` is
re-exported here. Do not import any of the following elsewhere in this
project -- the evaluation found them broken or miscalibrated:
  - edgeworth_upper_tail   (clips negative tail to 1e-300 instead of raising)
  - exact_median_pvalue    (out of scope for this pass; not adopted)
  - aggregate_network_null (55/66 analytic null means negative on the
                             package's own benchmark -- algebraically broken)
  - network_null_moments   (deferred to a later pass, not this one)

If HetNetEX-MD's own code ever needs a fix (e.g. one of the above), that is
a PR against lagillenwater/HetNetEx-MD, never a direct edit under
HetNetEX-MD/ in this repo.
"""

from __future__ import annotations

from hetnetex_md.core import exact_resampling_moments

__all__ = ["exact_resampling_moments"]
