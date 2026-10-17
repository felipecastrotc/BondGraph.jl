# 23: Causality assignment (SCAP) as a diagnostic

**What to build:** Assign causal strokes with the sequential causality assignment procedure (Karnopp §2.4, §3.4, Tables 3.6–3.7; Breedveld §1.4) and report what it finds: storage elements in derivative causality, algebraic loops, causal conflicts. Optionally draw the strokes in `bgplot`. This is the textbook analysis that plain acausal modelling doesn't give, and it would have warned about the pendulum's inertias in derivative causality before simulating (ticket 21).

Already decided (spec, "Out of scope"; ticket 10 of the map): causality is left to MTK; SCAP is an optional diagnostic and never replaces MTK's equations by default. BondGraphs.jl has the same item open (its issue #90).

**Blocked by:** None.

**Status:** needs-triage

Registered 2026-10-17 from the assessment "does the library add value?" (user asked to record the gaps and reassess later). Reassess before starting: does it still make sense, given MTK's own features and BondGraphs.jl's roadmap (`.scratch/mtk-upgrade/research/01-bondgraphs-jl-comparison.md`)?

## Comments
