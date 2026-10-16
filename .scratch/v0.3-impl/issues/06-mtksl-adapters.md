# 06: Standard-library adapters

**What to build:** A modeller attaches an MTK standard-library network to a bond through an adapter. The adapters are `ElectricalAdapter`, `TranslationalAdapter`, `TranslationalPositionAdapter`, `RotationalAdapter` and `HydraulicAdapter`. Each is a one-port with its law in bond variables, conventionally drawn into the adapter.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`. Prior art: research 11 (`../../mtk-upgrade/research/11-coupling-with-mtk-stdlib.md`) and its scripts.

**Blocked by:** 04.

**Status:** ready-for-agent

- [ ] ModelingToolkitStandardLibrary is a regular dependency with a `[compat]` bound, unless 05 already added it.
- [ ] The five adapters use the MTKSL port names (`p`/`n`, `flange`, `port`). The position-based ones take `s = nothing` / `phi = nothing`.
- [ ] The adapter signs are re-verified under the v0.2 convention (research 11 used the withdrawn convention A), and the docstrings state the sign rule and gotchas (the hydraulic fluid is connected in its own statement).
- [ ] Each adapter matches a pure-MTKSL reference model: RC, DC motor (through `RotationalAdapter` and through `ElectricalAdapter`), translational velocity and position, and hydraulic.
- [ ] A docs example per domain, or one page covering all five.
