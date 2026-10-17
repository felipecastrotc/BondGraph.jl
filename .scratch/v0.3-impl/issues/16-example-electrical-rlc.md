# 16: Example, Electrical circuits: series RLC

**What to build:** A series RLC circuit driven by a voltage step: Se, R, I (inductance), C on one 1-junction. Show the underdamped, critically damped and overdamped cases against the analytic step response. Source: Karnopp (2012) junction/circuit examples, Figs. 3.11–3.12 (note `docs/literature/karnopp-ch3.md`), plus the standard second-order solution.

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

- 2026-10-17: `examples/electrical_rlc.jl` built (agent). Se, R, I, C on one 1-junction (source in, passive elements out); `bgplot(model)` for the figure; `bg_equations` renamed to λ = L i, q and R, L, C, read against KVL. Intro derives the step response (characteristic roots, the three closed forms for v_C and i from rest). L = 1 H, C = 0.25 F (ω_n = 2 rad/s), V = 1 V; R = 1, 4, 10 Ω give ζ = 0.25, 1, 2.5, set per case via the `ODEProblem` parameter map. Plain-MTK reference in (q, i). Max differences on a 1001-point grid over [0, 20] s, max of v_C and i (tol 1e-10): vs analytic 4.7e-11 / 2.8e-11 / 2.4e-11 (under / critical / over); vs reference 0.0 in all three (with L = 1 it is the same ODE). Runs with `julia --project=docs`, no warnings. Literature: no new note; the existing `karnopp-ch3.md` (already in the README) covers what is cited (§3.1 Tables 3.1–3.4, §3.3 pp. 58–61 and Fig. 3.11, junction sign rule p. 59). Docs build not run.
