# How do power ports couple with ModelingToolkitStandardLibrary components?

Type: research
Status: open

## Question

With the accepted design (ADR 0001: `Power` connector with effort as the shared variable and flow as `[connect = Flow]`), how can a bond-graph model be coupled to ModelingToolkitStandardLibrary components: electrical `Pin` (v, i), translational `Flange` (s, f) / `MechanicalPort` (v, f), rotational `Flange` (φ, τ), hydraulic `HydraulicPort` (p, dm), and causal `Blocks` (`RealInput`/`RealOutput`) for sources and modulation? For each domain: is a direct `connect` possible (variable order/naming, which variable is the effort vs the flow, position- vs velocity-based connectors), or is a small adapter component needed, and what would it look like? Verify on MTK 11 with a runnable example per domain (e.g. a bond-graph DC motor driving an MTKSL rotational inertia; a bond-graph RC fed by an MTKSL voltage source; a `RealInput`-modulated GY).
