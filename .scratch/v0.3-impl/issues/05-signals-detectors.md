# 05: Signals: driven sources, signal moduli, detectors

**What to build:** A modeller drives a source or a two-port's modulus from MTK standard-library blocks (MSe, MSf, MTF, MGY). They measure an effort or a flow with a detector `De`/`Df`, which carries no power, and can close a feedback loop entirely in MTK.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 04.

**Status:** done (feature/mtk-upgrade)

- [x] ModelingToolkitStandardLibrary is a regular dependency with a `[compat]` bound, unless 06 already added it.
- [x] `Se()`/`Sf()` with no value expose `input::RealInput`. `TF()`/`GY()` with no modulus take it from `input::RealInput`.
- [x] `De` (zero flow, outputs the effort) and `Df` (zero effort, outputs the flow) expose `output::RealOutput`.
- [x] Check: each `input` form matches its value form driven by the same signal.
- [x] Check: a closed loop `De → Gain → Add → GY.input` (modulus `r = 2 + 0.5 e_C` read by a `De` on the 0-junction at C) matches the closed form `q = -4 + 2√(4 + 2t)`. (reworded 2026-10-17: the test uses `De`, see ticket 10)
- [x] A docs example for driven sources, signal moduli and detectors.

## Comments

- 2026-10-15: done. `@oneport`/`@twoport` add a signal connector the body defines as `input`/`output`. A modulus read from a signal leaves an algebraic equation (MTK won't divide by an unknown), so MTF/MGY driven by blocks need a DAE-capable solver (`Rodas5P`), plus `OrdinaryDiffEqNonlinearSolve` for the initial values of a nonlinear loop; both are test deps and the docs example (`examples/signals.jl`) explains it. The docs env needs MTKSL, `OrdinaryDiffEqRosenbrock` and `OrdinaryDiffEqNonlinearSolve` (ticket 10).
