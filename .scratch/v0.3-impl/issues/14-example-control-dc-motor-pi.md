# 14: Example, Control systems: DC motor speed control with a PI controller

**What to build:** A closed loop around the DC motor: `Df` measures the shaft speed, an MTKSL `Blocks` PI controller compares it with a reference and drives a modulated `Se` (the armature voltage). Shows signals in a feedback loop. Source: Michigan CTMS, *DC Motor Speed: PI(D) Controller Design* (same site and parameters as the DC motor example's link); compare the closed-loop response with the hand-derived closed-loop equations (plant + PI).

**Blocked by:** None (can start immediately). Follows the template of ticket 13.

**Status:** done (feature/mtk-upgrade), user review pending

## Shared decisions (user, 2026-10-17)

- One canonical example per application listed in `docs/src/index.md`, plus the hydraulic domain and a reaction network. Standard models that the literature uses as reference examples.
- Full ticket-13 template: physics and hand-derived equations with their source, bond graph and element ↔ parameter mapping, construction with each bond's direction, `bg_equations` read against the hand derivation, simulation, plain-MTK reference model, numeric comparison, conclusion.
- The bond graph figure is the library's own drawing (`bgplot(model)`), not a new TikZ file. No new image files; a physical sketch only as a link to the source.
- Source first: a short note in `docs/literature/` (paraphrase, what the source states vs **Inference:**, page/equation references) written only from what was actually read: the PDFs in `literature/` or an open-access source. No page or equation number that wasn't read.
- Textbook bond-graph conventions; never reinterpret the modeller.

## Acceptance criteria

- [x] Literature note written and linked from `docs/literature/README.md`.
- [x] Example follows the template; the comparison's max difference recorded below.
- [x] Runs with `julia --project=docs` without errors or initialisation warnings; the docs build.

## Comments

- 2026-10-17: Built `examples/control_dc_motor_pi.jl` and `docs/literature/ctms-dc-motor-speed.md` (agent). The live CTMS site answers with a bot check (HTTP 403), so both pages were read as Internet Archive snapshots (2022-11-19, 2023-06-10); figures and plot annotations not read. Loop: `Df` on the shaft junction → `Blocks.Feedback` (reference from `Blocks.Constant`, 1 rad/s) → `Blocks.PI` → modulated `Se()`; plant as the DC motor example without the torque source. CTMS's parameters (J 0.01, b 0.1, K 0.01, R 1, L 0.5) turn out to be the same as the DC motor example's; only the inputs differ (stated in the text). CTMS ends with a PID (Kp 100, Ki 200, Kd 10) and gives no PI; the example's PI gains are its own (Kp 20, Ki 40, zero on the slow plant pole), marked as such. Against the hand-derived closed loop as a plain-MTK model, 1001-point grid over 0–3 s: max difference 2.1e-10 (ω), 2.5e-10 (i), 4.2e-9 (V). Against CTMS's requirements: overshoot 1.7 % (< 5 %), 2 % settling 0.59 s (< 2 s), steady-state error 3e-4 % (< 1 %). Runs with `julia --project=docs`, no warnings; docs build not run; README link left to the user (first criterion open).
- Library gap (src unchanged): `bgplot` draws bonds only. The signal path Df → blocks → MSe isn't drawn, and the modulated source is labelled `Se`, not `MSe`. The example says so under the figure.
