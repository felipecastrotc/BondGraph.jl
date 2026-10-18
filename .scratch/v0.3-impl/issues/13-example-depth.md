# 13: Examples at the DC motor's level of detail

**What to build:** Every example explains its system the way the DC motor example does: the physics and the hand-derived equations, the bond graph and what each element stands for, a step-by-step construction, the equations the library generates set against the hand-derived ones, and a simulation compared with a plain-MTK model of the hand-derived equations. The v0.3 port cut the mechanical examples to about a third of their v0.2 length (user, 2026-10-17: "use it as template and add the same details", skipping a part only where it doesn't make sense for the example).

**Blocked by:** None (can start immediately).

**Status:** done (feature/mtk-upgrade), user review pending

## Template (from `examples/electrical_dc_motor.jl`)

1. Introduction: the physical system, its figure (circuit or free-body diagram), assumptions, the governing equations derived by hand (Newton, Kirchhoff), a reference.
2. Bond graph model: the TikZ figure, and which element stands for which physical parameter.
3. Preliminary steps: what the example does, the imports.
4. Element definitions, with units.
5. Model construction, with each bond's direction explained.
6. Analysing the model: `mtkcompile`, `unknowns`, `parameters`, the equations, and `bg_equations` read against the hand derivation (new, the DC motor gets it too).
7. Simulation of the bond graph model.
8. Reference model: the hand-derived equations as a plain MTK `System`, simulated.
9. Comparison plots and a one-line conclusion.

## Decisions

- Physical-model examples get the whole template: DC motor (adds step 6's `bg_equations`), harmonic oscillator, 2-DoF, N-DoF, Duffing.
- The Quick Example (RC) stays quick: it is spliced into Getting Started. It gets the hand-derived comparison only if it stays short.
- API guides (`two_ports`, `signals`, `adapters`) get the parts that fit: a hand-derived or closed-form reference and a comparison where they simulate a physical model; no physics introduction for an API feature.
- Out: `drawing`, `equations_latex`, `elements_oneport`. They are API pages without a system to derive by hand (`elements_oneport` already compares with the built-in element).
- Figures are not changed (the DC motor and 2-DoF TikZ fixes are separate). Text that disagrees with a figure is flagged, not silently worked around.

## Acceptance criteria

- [x] Each in-scope example follows the template, or says nothing about a skipped part that doesn't apply.
- [x] Each comparison is checked numerically when written (max difference recorded in a comment below), and the plots overlap.
- [x] Hand derivations follow the textbook sign conventions; what a book states is cited, the rest is derivation.
- [x] Every example runs with `julia --project=docs`; the docs build.

## Comments

- 2026-10-17: DC motor, RC, two-ports done (agent). DC motor gains `bg_equations` renamed to λ = L i, h = J θ̇ and the intro's symbols, read against the hand derivation; the intro now includes the applied torque τₐ the model has; citations for the GY law (Karnopp Eq. 3.11) and orientations (§3.1, Fig. 2.3d). RC: equation now in q (was V with charge in the text), closed-form overlay. Two-ports: hand derivation and closed form for the three lever runs. Max differences on 1001-point grids: DC motor 1.1e-10 (i), 5.4e-10 (ω); RC 6.7e-12; TF constant 3.6e-14, modulated 2.5e-10, custom 3.6e-14.
- Figure mismatches found (figures unchanged, for the TikZ fix): `bg_rc.tex` draws R→0 and C→0, code and text bond 0→R, 0→C; `bg_dcmotor.tex` draws jm→τ (known), labels the gyrator MGY though it is constant, and uses L_e, R_e, J_θ where the text uses L, R, J.
- 2026-10-17: 1-DoF, 2-DoF, Duffing done (agent): full template, free-body figures embedded, `bg_equations` renamed and latexified against the hand derivation, plain-MTK references. 1-DoF: v0.2's `c d²x/dt²` fixed to `c ẋ`. 2-DoF: the `for` loop (nothing rendered) replaced by a `two_dof(force; name)` builder and top-level calls. Max differences (unforced / forced): 1-DoF 0.0 / 1.1e-15; 2-DoF 1.4e-11 / 3.8e-11; Duffing 1.6e-10 / 6.5e-10.
- Figure mismatch (for the TikZ fix): `bg_2dof` draws the 0-junction bonds as forced mass → 0 → coupling and grounded mass; the code (and the regression data) has grounded mass → 0 → forced mass and coupling, so the coupling q is x₁ − x₂ in the code, x₂ − x₁ in the figure. Same physics; the text follows the code.
- 2026-10-17: N-DoF, signals, adapters done (agent). N-DoF: full template, Newton's law for the chain in general form and for N = 4; the 0-junction flow balance gives qᵢ = xᵢ₋₁ − xᵢ (compression positive), which is why the unit initial stretch is q₂ = −1 (the v0.2 remark is gone); no figure (text describes the bond graph). Signals: closed forms for the sine-driven RC (8.1e-10) and the feedback loop q = √(16+8t) − 4 (7.7e-13). Adapters: motor vs pure-MTKSL motor with `EMF` (ω 2.1e-11, i 1.7e-10), circuit vs pure-MTKSL (0.0 on the grid; 1.7e-11 vs 3(1−e⁻ᵗ)), translational vs v = t/2 (1.0e-14), hydraulic vs closed form (4.3e-12 relative; MTKSL's hydraulic resistances aren't linear, so no MTKSL reference). N-DoF max difference 1.9e-11.
- To check by the user: the 2-DoF text reads the coupling q = x₁ − x₂ as elongation (vertical stack, x down), the N-DoF text reads the same pattern's q as compression (horizontal chain). Both may be right for their own axes, but N-DoF points to 2-DoF as "the same pattern"; check against `fbd_2dof` and make the two consistent.
