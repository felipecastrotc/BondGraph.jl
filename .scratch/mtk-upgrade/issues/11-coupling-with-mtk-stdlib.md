# How do power ports couple with ModelingToolkitStandardLibrary components?

Type: research
Status: resolved

## Question

With the accepted design (ADR 0001: `Power` connector with effort as the shared variable and flow as `[connect = Flow]`), how can a bond-graph model be coupled to ModelingToolkitStandardLibrary components: electrical `Pin` (v, i), translational `Flange` (s, f) / `MechanicalPort` (v, f), rotational `Flange` (φ, τ), hydraulic `HydraulicPort` (p, dm), and causal `Blocks` (`RealInput`/`RealOutput`) for sources and modulation? For each domain: is a direct `connect` possible (variable order/naming, which variable is the effort vs the flow, position- vs velocity-based connectors), or is a small adapter component needed, and what would it look like? Verify on MTK 11 with a runnable example per domain (e.g. a bond-graph DC motor driving an MTKSL rotational inertia; a bond-graph RC fed by an MTKSL voltage source; a `RealInput`-modulated GY).

## Answer

See [research/11-coupling-with-mtk-stdlib.md](../research/11-coupling-with-mtk-stdlib.md) (ran on MTK 11.45.4 + MTKSL 2.29.8; scripts in its appendices).

- No MTKSL connector can be `connect`ed directly to `Power`: MTK errors on mismatched connector variable names (no silent mis-pairing). Renaming `Power`'s variables would not help (MTKSL mechanics' flow variable is force/torque = bond-graph *effort*; `Pin.v` is an absolute potential).
- A 2–3-equation **adapter element** works in every domain: the MTKSL side acts as the element at the end of the bond; through-variable negated; across-variable → effort (electrical, hydraulic) or, via its derivative, flow (mechanical). In a `Junction1` the adapter takes the sign of the element it replaces (+1 like a source, −1 like R/C/I).
- Verified against reference models: electrical RC (≤ 3.5e-10), DC motor ↔ MTKSL rotational (≤ 5.7e-10), translational velocity/position (≤ 7.7e-11), `RealInput`-modulated GY open/closed loop (bitwise / 7.7e-13), hydraulic (≤ 4e-11 rel).
- Recommendation (for the API ticket): five named adapters (`ElectricalAdapter`, `TranslationalAdapter`, `TranslationalPositionAdapter`, `RotationalAdapter`, `HydraulicAdapter`); `RealInput` forms of `Se`/`Sf`/`GY`/`TF`; `EffortSensor`/`FlowSensor`; MTKSL as a regular dependency (MIT, MTKBase-only, ~0.2 s load); consider renaming `Power` → `PowerPort` (clashes with `Blocks.Power`; `Mass`/`Spring`/`Damper` clash with MTKSL mechanics).
- Gotchas: `GY(g = blk.output.u)` needs `@named` (ParentScope); position adapters need an initial position; hydraulic needs the fluid `connect` and tiny `abstol`. Open: adapter sign rule needs maintainer sign-off; Thermal/Magnetic not covered.
