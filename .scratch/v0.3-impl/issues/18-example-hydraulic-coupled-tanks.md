# 18: Example, Chemical process / hydraulic: two coupled tanks

**What to build:** Two tanks in series with nonlinear outlet orifices (Torricelli, `Q = a√(2 g h)` or in pressure form `Q ∝ √p`), fed by a pump (flow source). Tanks are C elements (pressure = ρ g V / A), orifices a nonlinear R written as a custom `@oneport`. The standard liquid-level process plant. Source: Borutzky 2011, Ch. 4 §4.4 coupled tanks (read the PDF, `literature/Wolfgang_Borutzky.pdf`); compare with the hand-derived level ODEs, and the steady-state levels in closed form.

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

- 2026-10-17 (agent): `examples/hydraulic_coupled_tanks.jl` built, following Borutzky Fig. 4.8 (p. 146): Sf:pump → 0 (C:tank1) → 1 (valve1) → 0 (C:tank2) → 1 (valve2) → Se:atm (0 Pa gauge), every bond downstream, Se:0's bond into it as drawn. Tanks are linear `Compliance(C = A/(ρg))`, volume as state; valves a custom `@oneport Orifice` with Eq. 4.12 in its sign-safe form `cd a sign(e) √(2|e|/ρ)` (p. 141). Nonlinear valves are our modelling choice: the book linearises them (Eq. 4.28). Tank 2 starts fuller, so valve 1's flow reverses (Q₁ < 0 until t ≈ 27 s), which exercises the sign-safe form. Literature: added a "Coupled hydraulic tanks in detail" section to `docs/literature/borutzky-ch3-4.md` (already linked from the README; README not edited). The book gives no parameter values for this example; ours are A = 0.1 m², a = 5 cm², cd = 0.6, Qp = 1 L/s. Max differences against the plain-MTK level ODEs (Torricelli form, no ρ), 1001-point grid over 0–5000 s, reltol 1e-10 / abstol 1e-12: h₁ 4.5e-9 m, h₂ 4.4e-9 m. Steady state vs closed form (h₁* = 1.1326 m, h₂* = 0.5663 m) at t = 5000 s: −2.3e-8 m, −1.4e-8 m. Runs with `julia --project=docs`, no errors or initialisation warnings; docs build not run. Library note: `bgplot` labels custom elements by their constructor name (`Orifice:valve1`, not the book's `R:R1`); only the member order sets the drawing, so `jv2` lists `atm` first to put it on the path, as in the book.
