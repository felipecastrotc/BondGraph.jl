# Prototype the bond connection mechanism on MTK 11

Type: prototype
Status: resolved
Blocked by: 03

## Question

Ticket 03 recommends the hybrid design (0-junction = native `connect` on a `Power` port with `Flow` flow; 1-junction = component with ±1-signed ports; TF/GY as two-port components) and already ran RC, DC motor and 2-DoF against the golden data. Does it hold for the remaining examples (N-DoF, Duffing, a modulated two-port, a custom nonlinear `@oneport`-style element), and does model-building code written with it keep today's ergonomics? A throwaway prototype (starting from the scripts referenced in research/03) should express all six examples in the proposed user-facing API and reproduce `test/regression/reference/`, so the maintainer can react to how model-building code looks with it.

## Answer

**Verdict (maintainer, 2026-10-15): accept the hybrid design.** Recorded as [ADR 0001](../../../docs/adr/0001-junctions-on-native-mtk-connections.md).

Prototype: branch `prototype/connection-mechanism` (local, never merged), `prototype/connection_mechanism/{BGProto.jl,run.jl,README.md}`. On MTK 11.45.4 all six examples (9 golden cases) reproduce `test/regression/reference/` (max |Δ| ≤ 2.3e-10); modulated GY conserves power exactly.

Decided with it:
1. 0-junction = native `connect`; 1-junction = signed-port core component; TF/GY = two-port components; no MTK internals.
2. Junction-to-junction bonds use named free ports (`:out`, `[-1, :out]`) and an explicit `connect(a.out, b.in, c.in)` per 0-junction (replaces v0.2's shared `.power` hub + argument-order direction).
3. Inertia keeps momentum `p` as its state with an initial-value keyword (`Mass(m, p = 0.0)`, mirroring `Spring(k, x)`); flow stays available as `power.f`.
4. `mGY(...; coneqs)` and the `@named sys, con = ...` tuple form are replaced by `GY(g = ...)` / `TF(n = ...)` components wired with `connect`; modulated = pass an expression.

Findings for other tickets:
- **Sign convention (06):** the N-DoF v0.2 model bonds the spring-damper *into* the 0-junction; bonding it out gives the mirrored physics (velocities negated). Bond direction choices are visible and matter.
- The v0.2 N-DoF "non-adjacent merge" anomaly does not occur in the new mechanism. Duffing's redundant state (`s.q`, `s3.q`) remains.
- MTK 11 discovers variables itself: `@oneport` shrinks to a 5-line convenience; `isindependent`-style sorting is unnecessary.
