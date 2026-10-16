# 05: Signals: driven sources, signal moduli, detectors

**What to build:** A modeller drives a source or a two-port's modulus from MTK standard-library blocks (MSe, MSf, MTF, MGY). They measure an effort or a flow with a detector `De`/`Df`, which carries no power, and can close a feedback loop entirely in MTK.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 04.

**Status:** ready-for-agent

- [ ] ModelingToolkitStandardLibrary is a regular dependency with a `[compat]` bound, unless 06 already added it.
- [ ] `Se()`/`Sf()` with no value expose `input::RealInput`. `TF()`/`GY()` with no modulus take it from `input::RealInput`.
- [ ] `De` (zero flow, outputs the effort) and `Df` (zero effort, outputs the flow) expose `output::RealOutput`.
- [ ] Check: each `input` form matches its value form driven by the same signal.
- [ ] Check: a closed loop `Df → Gain → Add → GY.input` matches `GY(r = expression)`.
- [ ] A docs example for driven sources, signal moduli and detectors.
