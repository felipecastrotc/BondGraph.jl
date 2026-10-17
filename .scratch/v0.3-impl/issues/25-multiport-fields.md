# 25: Multiport fields (C, I, R with several ports)

**What to build:** Multiport storage and dissipation fields: C- and I-fields defined by an energy function, R-fields, e.g. a two-port C with `E = q₁²/(2 C(q₂))` (Breedveld §1.5.11–1.5.12, Eqs. 1.44–1.49; Karnopp's later chapters, not yet read). Today only one-ports (`@oneport`) and two-ports (`@twoport`) exist; a field can be written by hand as an MTK component, but there is no helper and no docs. Needed for electromechanical transducers and coupled storage.

**Blocked by:** None.

**Status:** needs-triage

Registered 2026-10-17 from the assessment "does the library add value?" (user asked to record the gaps and reassess later). Reassess before starting: does it still make sense, given MTK's own features and BondGraphs.jl's roadmap (`.scratch/mtk-upgrade/research/01-bondgraphs-jl-comparison.md`)?

## Comments
