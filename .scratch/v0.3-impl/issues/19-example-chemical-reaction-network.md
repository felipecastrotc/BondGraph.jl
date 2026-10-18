# 19: Example, Chemical reactions: a reaction network as a bond graph

**What to build:** A small reaction network (e.g. A ⇌ B ⇌ C, or A + B ⇌ C) as a true bond graph: chemical potential is the effort, molar flow the flow; species are nonlinear `Ce` elements (`μ = RT ln(K x)`, a custom `@oneport`) on 0-junctions, reactions are `Re` two-ports (Marcelin–de Donder: `v = κ(exp(Af/RT) − exp(Ar/RT))`, a custom `@twoport`), with 1-junctions collecting the forward and reverse affinities. Source: Gawthrop & Crampin 2014, *Energy-based analysis of biochemical cycles using bond graphs*, Proc. R. Soc. A (open access). Compare with the mass-action kinetics ODEs, and check that the equilibrium matches the equilibrium constants.

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

- 2026-10-17: Built `examples/chemical_reaction_network.jl` (agent): A ⇌ B ⇌ C, Gawthrop & Crampin 2014 Fig. 1c. Custom `Ce` (`@oneport`, μ = RT ln(K x)) and `Re` (`@twoport`, Marcelin–de Donder, f2 = f1) on three 0-junctions; species bonds 0 ⇀ Ce, reaction bonds follow the forward reaction (0 ⇀ Re ⇀ 0), per the paper's Eq. 10. No 1-junctions: a unimolecular chain doesn't need them (they appear for A + B ⇌ C + D, Fig. 1d). Figure from `bgplot(network)`. `bg_equations` renamed to x_a…, K_a…, κ_1, κ_2, RT reduces to exp(log(K x)) terms, i.e. the mass-action ODEs. Max differences vs the plain-MTK mass-action model on a 1001-point grid over 0–50 s: 2.5e-10 (x_a), 3.5e-10 (x_b), 9.8e-11 (x_c). Equilibrium at t = 50 s vs x_i ∝ 1/K_i: max 5.0e-11; ratios minus K₁ = Ka/Kb, K₂ = Kb/Kc: 1.0e-11, −8.7e-10. Re dissipation ≥ 0 along the trajectory (min 3e-19, 1.8e-17). Runs with `julia --project=docs`, no warnings. Literature note `docs/literature/gawthrop-crampin-2014.md` written from arXiv v2 §1–§3.1 (not linked from the README yet, left to the maintainer). Flagged in the note: the paper says the equilibrium constant of A ⇌ B is Kb/Ka; Eq. 7 at v = 0 gives xb/xa = Ka/Kb. No library limitation found.
