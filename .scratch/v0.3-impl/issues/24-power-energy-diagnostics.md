# 24: Power and energy diagnostics

**What to build:** Power and energy checks from the bond variables every bond already has: power on a bond (`e·f`), the power balance at each junction and two-port (Karnopp Eq. 3.16), stored energy of C and I elements and dissipation in R, and the activity index used for model reduction (Louca et al. in Borutzky Ch. 2, Eqs. 2.2–2.6; note `docs/literature/borutzky-ch2-reduction.md`). Useful for verifying models and for teaching. Cheap to build. The reaction-network example already checks reaction power by hand.

**Blocked by:** None.

**Status:** needs-info

Registered 2026-10-17 from the assessment "does the library add value?" (user asked to record the gaps and reassess later). Reassess before starting: does it still make sense, given MTK's own features and BondGraphs.jl's roadmap (`.scratch/mtk-upgrade/research/01-bondgraphs-jl-comparison.md`)?

## Comments

- 2026-10-18 (triage, user agreed): after v0.3. Reassess against MTK's features and BondGraphs.jl's roadmap before specifying.
  Cheapest of the three value gaps: first in line for 0.4.
