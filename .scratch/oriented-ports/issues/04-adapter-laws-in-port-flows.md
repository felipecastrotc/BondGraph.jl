# 04: An adapter bonded through `pout` breaks power conservation

**What to build:** an adapter is a power-conserving coupling to an MTKSL network, not a constitutive one-port. Its laws are written in the bond's `f` ("current into the network"), so under ADR 0002's "you get what you drew" an adapter bonded through `pout` counts the network's power with the wrong sign on the bond-graph side. Proposed fix: write the network-side law in the ports' MTK flows, e.g. `ElectricalAdapter`'s current into the network `= power.f + pout.f`, so either port conserves power and `adapter.pout` means "the network supplies power".

**Blocked by:** None.

**Status:** done (feature/mtk-upgrade, uncommitted), user review pending

Found by rechecking usability trial 02 (`.scratch/api-usability/runs/02/`, finding P02-F01 / U03) on the new API.

## Reproduction

MTKSL source (1 V, R = 2, L = 0.2) supplying a BG motor (`GY(r = 0.1)`, I = 0.02, R = 0.02 + 0.03, C = 2), checked against Kirchhoff/Newton:

| Interface | Result |
|---|---|
| `Junction0(adapter, :motor)`, `connect(interface.motor, motor.p1)` (the trial's) | `Success`, max error 1.3e-10 |
| `connect(adapter.pout, motor.p1)` | `Success`, drawn `adapter ⇀ GY ⇀ 1`, **max error 0.64** |

No error or warning: the natural form gives a wrong result. Script: `.scratch/api-usability/recheck/trials_02_05.jl` (`julia --project=docs`).

## Decision needed

- Does ADR 0002's "laws in `e`, `f`" apply to adapters? The fix exempts them: an adapter's MTKSL side is a network, not a law, so the direction of its bond can't change the physics. Amend ADR 0002 (consequences) if accepted.
- Alternative: keep the laws and document that adapters take only `power` (reject `pout` in `bondgraph`, like the both-ports error). Weaker: still silent at `mtkcompile`.

## Acceptance criteria

- [x] A supplying MTKSL network through `adapter.pout` matches the Kirchhoff/Newton reference (trial 02 numbers), for all five adapters.
- [x] The receiving direction (`power`) results are unchanged; existing tests and examples pass.
- [x] `examples/adapters.jl` shows the supplying direction (closes U03).

## Comments
- 2026-10-18 (maintainer): adapters work in both directions (port-flow laws).
- 2026-10-18 (agent): network-side laws now use `power.f + pout.f` (`src/adapters.jl`); effort laws unchanged. Tests in `test_adapters.jl`: the trial-02 armature through `pout` and through the old 0-junction form, electrical and hydraulic sources through `adapter => :pout` (positive charge/volume, current out of the network), and an MTKSL force/torque on a bond-graph inertia through each mechanical adapter's `pout`. Full suite 174/174; `examples/adapters.jl` runs, its hydraulic tank now bonds through `pout` (volume positive, closed form error 4e-12 relative). Recheck script: `pout` error 1.4e-10 (was 0.64). ADR 0002 amended. Consequence to know: in mechanical domains the effort stays the force on the network, so through `pout` the bond-graph side's velocity is the flange's mirrored (pinned by the test, documented in the example table).
