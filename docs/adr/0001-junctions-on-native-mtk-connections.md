# 1. Junctions on native ModelingToolkit connections

Date: 2026-10-15
Status: accepted

## Context

BondGraph v0.2 (MTK 8.71) turned bonds into equations with its own algorithm: it overrode MTK internals (`expand_connections`, `generate_connection_set`, `@named` tuple handling) to pull bond-graph connections out of MTK's connection sets, kept `connect` argument order as bond direction, and generated junction and two-port equations itself. MTK 10.4 rewrote connection expansion, and Julia 1.10+ refuses method overwriting during precompilation, so that design cannot be ported. Research showed MTK 11's native acausal `connect` already has 0-junction semantics.

## Decision

- **0-junction = MTK's native `connect`** on a single `Power` connector (effort shared, flow `[connect = Flow]`, positive into the component).
- **1-junction = a component** with one signed `Power` port per bond (+1 into, −1 out of the junction).
- **TF / GY (and modulated forms) = ordinary two-port components**, wired with `connect`; modulation is just an expression for the ratio.
- No overrides of MTK internals; simulation is plain `mtkcompile` → `ODEProblem` → `solve`.

## Consequences

- The custom connection/equation-generation machinery (`src/mtk.jl`, most of `src/utils/`, bond metadata in `src/types.jl`) is deleted, not ported.
- Bond direction is no longer `connect` argument order: element bonds point into the element; junction-to-junction orientation is the ±1 on 1-junction ports.
- A 1-junction needs one port per bond (named free ports), not a shared hub port; `mGY(...; coneqs)` and the `@named sys, con = ...` tuple form go away.
- Evidence: a prototype (branch `prototype/connection-mechanism`) reproduced all six v0.2 examples' golden data on MTK 11.45.4 (max |Δ| ≤ 2.3e-10).
