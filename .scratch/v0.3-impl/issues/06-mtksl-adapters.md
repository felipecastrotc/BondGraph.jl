# 06: Standard-library adapters

**What to build:** A modeller attaches an MTK standard-library network to a bond through an adapter. The adapters are `ElectricalAdapter`, `TranslationalAdapter`, `TranslationalPositionAdapter`, `RotationalAdapter` and `HydraulicAdapter`. Each is a one-port with its law in bond variables, conventionally drawn into the adapter.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`. Prior art: research 11 (`../../mtk-upgrade/research/11-coupling-with-mtk-stdlib.md`) and its scripts.

**Blocked by:** 04.

**Status:** done (feature/mtk-upgrade)

- [x] ModelingToolkitStandardLibrary is a regular dependency with a `[compat]` bound, unless 05 already added it.
- [x] The five adapters use the MTKSL port names (`p`/`n`, `flange`, `port`). The position-based ones take `s = nothing` / `phi = nothing`.
- [x] The adapter signs are re-verified under the v0.2 convention (research 11 used the withdrawn convention A), and the docstrings state the sign rule and gotchas (the hydraulic fluid is connected in its own statement).
- [x] Each adapter matches a pure-MTKSL reference model: RC, DC motor (through `RotationalAdapter` and through `ElectricalAdapter`), translational velocity and position, and hydraulic.
- [x] A docs example per domain, or one page covering all five.

## Comments

- 2026-10-15: done. Research 11's adapter equations hold under the v0.2 convention with the adapter's bond always drawn into it (`[-1, adapter]`): MTKSL loads (electrical R+C, rotational inertia+damper behind the DC motor, translational masses) match their references exactly; MTKSL sources (voltage, pressure) give the reference with the opposite sign, i.e. negative power along the bond, as the convention says. Position adapters take `s`/`phi` like MTKSL components (a local variable with a `nothing` default). Docs example: `examples/adapters.jl`.
