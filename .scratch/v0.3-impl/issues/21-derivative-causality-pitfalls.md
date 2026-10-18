# 21: Derivative-causality models: initial values and MTK's state choice

**What to build:** Two pitfalls found by the ticket-15 example (pendulum driven by a DC motor), where two inertias are in derivative causality and MTK reduces the index:

1. The inertias in derivative causality need `Inertia(p = nothing)`, or initialisation is over-determined and warns. It works, but nothing documents it.
2. MTK picks `p_x` as a state and recovers ω from `p_x = −m l sin φ · ω`, which is singular at `sin φ = 0`. Started at φ = 0, Rodas5P reports success while the pendulum stays frozen: a silent wrong answer. The example avoids it by staying in −π < φ < 0 and says so.

Decide: document both (where?), and/or let the modeller steer MTK's state choice (e.g. MTK's `state_priority` metadata) in the example or in the elements.

**Blocked by:** None.

**Status:** done (feature/mtk-upgrade, uncommitted), user review pending

## Comments

- 2026-10-18 (triage, user agreed): **v0.3**, docs only: document `Inertia(p = nothing)` (and `Compliance(q = nothing)`) for elements in derivative causality, and the silent frozen pendulum at sin φ = 0 (docstrings + user guide). Try MTK's `state_priority` in the pendulum example; if it steers the state choice to φ, ω, use it there and say so. No src/ behaviour change.
- 2026-10-18 (agent): `state_priority` works. With the built-in inertias MTK keeps `p_x`; started at φ = 0 the run returns Success but is off by 1.33 rad (frozen). With `@variables p(t) [state_priority = -10]` on the point mass (a small `@oneport` in the example, since `Inertia` has no such keyword and src behaviour stays unchanged) MTK keeps the rotor momentum, nothing divides by sin φ, and a new run from φ = 0 matches the reduced model (7.5e-12). Docs: `Inertia`/`Compliance` docstrings (`p`/`q = nothing`), user guide "Simulation and causality", the example's analysis rewritten for the new states. Possible later: a `state_priority` keyword on `Inertia`/`Compliance` (fits ticket 27).
