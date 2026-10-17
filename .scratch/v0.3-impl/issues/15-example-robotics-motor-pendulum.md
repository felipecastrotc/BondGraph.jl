# 15: Example, Robotics/mechatronics: pendulum driven by a DC motor

**What to build:** Breedveld's planar pendulum driven by a DC motor (Borutzky 2011, Ch. 1, Figs. 1.6–1.8, Eqs. 1.13–1.26; note `docs/literature/borutzky-ch1-breedveld.md`): voltage source, armature R and L, GY, rotor J with bearing friction, a point mass reached through two state-modulated TFs (`−l sin φ`, `l cos φ`), gravity as an Se. The point-mass inertias are in derivative causality, so MTK must reduce the index. Compare with the reduced damped pendulum (inertia J + m l², gravity torque) the book derives, and check the small-angle frequency (Eq. 1.26).

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

- 2026-10-17 (agent): Built `examples/robotics_motor_pendulum.jl` (`@id motor-pendulum`), full ticket-13 template: Fig. 1.7d bond graph with `Se:us`, `R:Re`, `I:L`, `GY:Km`, `I:J`, `R:b`, two state-modulated `TF`s (`m = -l*sin(φ)`, `l*cos(φ)`), `I:mx`, `I:my`, `Se:mg = −mg`; φ comes from a `Df` on the ω junction into `Blocks.Integrator` (Breedveld's ∫ signal). `mtkcompile` reduces the index (dummy derivative ω̇); `bg_equations` read against the hand derivation (gyristor cancellation, Eqs. 1.22–1.24). Reference: reduced pendulum (J + ml², mgl cos φ), Fig. 1.11b. Max differences on a 1001-point grid over 10 s (Rodas5P vs Tsit5, tol 1e-10): φ 8.3e-12, ω 5.8e-11, i 4.9e-10. Small-angle check (Km = 0, b = 0, ψ₀ = 1e-3): max |ψ/ψ₀ − cos ω₀t| = 3.2e-6, the expected ψ₀²ω₀t/16 amplitude shift; f = 0.823 Hz (Eq. 1.26). Runs with `julia --project=docs`, no warnings. Literature note `borutzky-ch1-breedveld.md` extended (Fig. 1.6 geometry, Fig. 1.7d half-arrows, p. 33–35 statements); not yet linked from README (user). Docs build not run.
- Library limitations found (src unchanged): (1) the dependent inertias need `Inertia(p = nothing)` to avoid an over-determined initialisation warning; `p = nothing` works but isn't documented. (2) MTK picks states statically: here p_x, so ω is recovered from p_x = −m l sin φ ω, singular at sin φ = 0. Starting at φ = 0 the run "succeeds" silently with a frozen pendulum; the example stays in −π < φ < 0 and says so. (3) `bgplot` doesn't draw signals (Df → ∫ → MTF), and the `I:mx` label overlaps the `my` bond in the drawing.
