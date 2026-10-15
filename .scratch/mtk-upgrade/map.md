# Map: MTK 11 upgrade of BondGraph.jl

Labels: wayfinder:map
Status: open

## Destination

A settled, written spec for BondGraph v0.3 on ModelingToolkit 11 / Julia ≥ 1.10 (installable with fresh dependency resolution), in which bond-graph models remain ordinary MTK systems: they compose and couple with any MTK model, simulate with standard MTK functions, and keep custom nonlinear elements and modulated two-ports as easy to write as today. Every design decision is recorded, ready for `/to-spec` → `/to-tickets` → `/implement`.

## Notes

- **Domain:** bond-graph modeling as an MTK dialect. Read `GLOSSARY.md` before any ticket and use its terms; update it via `domain-modeling` when terms settle. ADRs go in `docs/adr/`.
- **Skills:** grilling tickets call `grilling` + `domain-modeling`; research tickets call `research`; prototype tickets call `prototype`.
- **Inputs:** `.scratch/mtk-upgrade/research-mtk-internals.md` (what changed in MTK 8.71 → 11 and why `src/mtk.jl` breaks), `.scratch/roadmap.md`, `TODO.md` (open sign-convention doubts), `.scratch/golden-tests/report.md` (regression safety net).
- **Invariant (from the maintainer):** keep the MTK-native ergonomics, i.e. tight integration and interchangeability with MTK models and functions, and easy custom ports, nonlinear R/C, and mGY/mTF. Everything else may break (no known users besides the maintainer; v0.2 → v0.3).
- **Physics must not change:** the six examples must reproduce `test/regression/reference/*.csv` (state names may be remapped in `test/regression/helpers.jl`).
- **Long-term goals that constrain decisions:** register the package in the General registry and publish a software paper (thesis contribution). The name must be registrable and the license paper-compatible.
- **Codex:** may take AFK work, but the maintainer approves every Codex launch.

## Decisions so far

<!-- one line per resolved ticket: [title](issues/NN-slug.md): gist -->

- [How does BondGraphs.jl compare today?](issues/01-bondgraphs-jl-comparison.md): complementary, not a replacement (graph object → MTK, biochemistry, no GY/modulated two-ports); `BondGraph` is not registrable, `BondGraphToolkit` and similar pass.

## Not yet specified

- **Coupling with non-bond-graph MTK components**, e.g. ModelingToolkitStandardLibrary electrical `Pin` / mechanical `Flange`: how a power port meets them, implied by the "couple with any MTK model" goal. Sharpens once the connection mechanism is chosen.
- **Causality**: how integral/differential causality and causal conflicts surface when MTK does index reduction; ties to the MTK vs ModelingToolkitBase stance.
- **Graph visualization**: rebuilding `generate_graph` on the new connection representation.
- **Regression oddities**: N-DoF 0-junctions merging non-adjacent elements (cancelling terms) and Duffing's redundant state; whether the new mechanism should fix or preserve them.
- **`Integral` / DomainSets** (`It` in `src/bg.jl`): keep, replace, or drop.

## Out of scope

- Documentation overhaul (roadmap step 4): follow-up effort once the API settles.
- Registration in the General registry: later, but the name decision must keep it possible.
- Software paper: later; informs positioning only.
- 2D / multibond experiments (`src/mbg.jl`): not in v0.3; delete during implementation.
