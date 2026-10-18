# 1. Junctions on native ModelingToolkit connections

Date: 2026-10-15
Status: accepted; amended 2026-10-15 by ticket 06 (sign convention) and 2026-10-18 by [ADR 0002](0002-oriented-ports.md) (element ports, junction member syntax)

## Context

BondGraph v0.2 (MTK 8.71) turned bonds into equations with its own algorithm: it overrode MTK internals (`expand_connections`, `generate_connection_set`, `@named` tuple handling) to pull bond-graph connections out of MTK's connection sets, kept `connect` argument order as bond direction, and generated junction and two-port equations itself. MTK 10.4 rewrote connection expansion, and Julia 1.10+ refuses method overwriting during precompilation, so that design cannot be ported. Research showed MTK 11's native acausal `connect` already has 0-junction semantics.

## Decision

- **0- and 1-junctions = one signed component** with one `PowerPort` per bond. The sign is the bond's half-arrow: +1 into the junction, −1 out of it. The equations are the textbook ones in bond variables (Karnopp Eqs. 3.16–3.21): the shared variable is unsigned and the sum is signed.
- **Element laws are written in the variables of the element's bond, as drawn** (v0.2's convention, unchanged). A junction converts MTK's port flow (positive into the component) to the bond flow, so a bond drawn out of an element keeps the element's law as written.
- **Plain `connect` of free junction ports is a 0-junction** whose bond directions come from the signs at the other end of each bond.
- **TF / GY (and modulated forms) = ordinary two-port components**, wired with `connect`; modulation is just an expression for the ratio.
- No overrides of MTK internals; simulation is plain `mtkcompile` → `ODEProblem` → `solve`.

## Consequences

- The custom connection/equation-generation machinery (`src/mtk.jl`, most of `src/utils/`, bond metadata in `src/types.jl`) is deleted, not ported.
- Bond direction is no longer `connect` argument order: it is the ±1 the modeler gives each junction argument, at both junction types. Junction-to-junction bonds carry their direction on the junctions' free ports.
- A 1-junction needs one port per bond (named free ports), not a shared hub port; `mGY(...; coneqs)` and the `@named sys, con = ...` tuple form go away.
- Evidence: a prototype (branch `prototype/connection-mechanism`) reproduced all six v0.2 examples' golden data on MTK 11.45.4 (max |Δ| ≤ 2.3e-10). After the ticket 06 amendment it still does, with v0.2's own signs. A resistor drawn reversed behaves as a negative resistor, and the textbook N-DoF orientation (`0 → sd`) matches the v0.2 data with the coupling springs' `q` negated.
- The first version of this ADR made the 0-junction plain `connect` and fixed every element's bond to point into the element. That ignored the signs modelers give at 0-junctions and reinterpreted reversed bonds. It was rejected in ticket 06: the tool follows bond-graph conventions and never assumes what the modeler meant.
