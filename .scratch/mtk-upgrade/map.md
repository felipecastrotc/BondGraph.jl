# Map: MTK 11 upgrade of BondGraphToolkit.jl

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
- **License deferred:** the dependency/license ticket is deferred by the maintainer; provisionally depend on full ModelingToolkit and keep the current license. Don't block other tickets on it.
- **Codex:** may take AFK work, but the maintainer approves every Codex launch.

## Decisions so far

<!-- one line per resolved ticket: [title](issues/NN-slug.md): gist -->

- [How does BondGraphs.jl compare today?](issues/01-bondgraphs-jl-comparison.md): complementary, not a replacement (graph object → MTK, biochemistry, no GY/modulated two-ports); `BondGraph` is not registrable, `BondGraphToolkit` and similar pass.
- [Can junctions be expressed with MTK's native connections?](issues/03-junctions-as-native-connections.md): yes, via a hybrid: 0-junction = native `connect` (effort Equality, flow `Flow`), 1-junction = component with ±1-signed ports, TF/GY = two-port components; no MTK internals; RC/DC motor/2-DoF match golden data on MTK 11.
- [Which MTK 11 dependencies are AGPL, and what does that mean for us?](issues/02-licensing-facts.md): full MTK ≥ 11 hard-depends on two AGPL-3.0 packages; ModelingToolkitBase is MIT but lacks index reduction/tearing (higher-index bond graphs fail); Catalyst's pattern (depend on MTKBase, users opt into full MTK) is the precedent; JOSS is fine either way.
- [Prototype the bond connection mechanism on MTK 11](issues/04-prototype-connection-mechanism.md): hybrid accepted ([ADR 0001](../../docs/adr/0001-junctions-on-native-mtk-connections.md)); all six examples reproduce golden data on MTK 11; free ports + explicit `connect`, momentum state for inertias, GY/TF components replace `mGY(coneqs)`.
- [Choose the connection mechanism](issues/05-connection-mechanism.md): settled by the ticket 04 verdict (ADR 0001).
- [How do power ports couple with ModelingToolkitStandardLibrary components?](issues/11-coupling-with-mtk-stdlib.md): no direct `connect`, but small per-domain adapter elements match reference models in electrical, translational, rotational, hydraulic, and `RealInput`-modulated cases; API choices go to the public-API ticket.
- [Fix the sign convention and the meaning of bond direction](issues/06-sign-convention.md): v0.2's convention stays, and the textbooks agree (half-arrow signs only in junction sums; one-port laws in bond variables; TF/GY port 1 in, port 2 out). 0- and 1-junctions become one signed component, and ADR 0001 is amended. The prototype passes all golden cases with v0.2 signs.
- [Choose the package name](issues/08-package-name.md): `BondGraphToolkit` (passes the registry checks, matches the existing docs). Module, `Project.toml`, file, docs and tests are renamed in v0.3 and the UUID is kept. The GitHub repo is renamed at the v0.3.0 release, after asking.
- [Decide the public API and its breaking changes](issues/10-public-api.md): role names `Resistance`/`Compliance`/`Inertia` (book parameters R, C, I; states q, p) plus `Mass`/`Spring`/`Damper` shorthands; `PowerPort`; `TF(m)`/`GY(r)` with expression or `RealInput` modulus; `De`/`Df`; `[-1, x]` junction syntax kept; `@oneport` + `@twoport`; five MTKSL adapters; `bg_equations` + `latexify`; Graphviz plot + TikZ export; causality left to MTK.
- [Confirm the target versions and compat bounds](issues/09-target-versions.md): MTK 11, `ModelingToolkit = "11.45"` (caret), Julia ≥ 1.10. Direct deps only for what the code imports. CI on Julia 1.10 and 1.

## Not yet specified

Nothing left: causality (left to MTK), graph visualization (Graphviz + TikZ), Duffing (kept as drawn) and `Integral`/DomainSets (dropped) were settled in ticket 10.

## Future (after v0.3)

- SCAP causality assignment as an optional diagnostic (reports derivative causality, conflicts and algebraic loops by element; strokes on the plot), never replacing MTK's equations by default.
- Multiport fields (`@multiport`), thermal/magnetic adapters.

## Out of scope

- Documentation overhaul (roadmap step 4): follow-up effort once the API settles.
- Registration in the General registry: later, but the name decision must keep it possible.
- Software paper: later; informs positioning only.
- 2D / multibond experiments (`src/mbg.jl`): not in v0.3; delete during implementation.
