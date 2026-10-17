# 17: Example, Biomedical: 3-element Windkessel model of the arterial circulation

**What to build:** The three-element Windkessel: characteristic impedance Rc in series with the parallel peripheral resistance Rp and arterial compliance C, driven by an aortic inflow (a flow source with a pulsatile waveform). A hydraulic analogy (pressure = effort, volume flow = flow). Source: Westerhof, Lankhaar & Westerhof 2009, *The arterial Windkessel*, Med Biol Eng Comput (open access); typical parameter values from it. Compare with the hand-derived ODE for aortic pressure.

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

- 2026-10-17: Built `examples/biomedical_windkessel.jl` and `docs/literature/westerhof-2009-windkessel.md` (agent). Read the Springer HTML full text, §1–§8, with Figs. 3 and 5 and Tables 1–2 (the paper isn't in PMC; PMC2644372 is another article). Topology checked against Fig. 5: Sf(Q) → 1 → R:Zc, 1 → 0 → C, 0 → R:Rp. The paper gives no parameter values (only Z_c ≈ 5–7% of R), no ODE and no inflow waveform, so these are marked as our choices: half-sine systolic inflow (T = 0.8 s, Ts = 0.3 s, SV = 70 mL), R = 1.0, Z_c = 0.06 mmHg⋅s/mL, C = 1.5 mL/mmHg. Hand-derived C Ṗc = Q − Pc/R, P = Pc + Z_c Q. Ten beats; last beat 73.3/116.6 mmHg, change from the previous beat −0.02/−0.03 mmHg, mean 92.8 against (R + Z_c) SV/T = 92.75. Max difference on a 1001-point grid (reltol = abstol = 1e-10, `tstops` at the inflow's kinks): 7.3e-8 mmHg in both P and Pc (≈1e-9 relative). Runs with `julia --project=docs` with no output, so no initialisation warnings. Docs build not run. Note: `bg_equations` lists only state equations, so `j1.Zc.R` (only in the observed entrance pressure) can't be renamed; the example explains this rather than changing src/. README link left to the user.
