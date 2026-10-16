# 07: Readable equations: `bg_equations` and `latexify`

**What to build:** A modeller gets the model's equations in terms of its states, renames variables for display (including inside derivatives, e.g. `D(s.q)` → `ẋ`, `m.p` → `p`), and exports them to LaTeX for a report or paper.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** done (feature/mtk-upgrade)

- [x] `bg_equations(sys; rename = Dict())` accepts a compiled or uncompiled model and returns MTK's full equations in terms of the states. The renames are applied everywhere, including inside `Differential`. The model is not modified.
- [x] `latexify` is re-exported from Latexify, which is a direct dependency with a bound.
- [x] Checks: renamed variables appear and the original names don't, including in derivatives; `latexify` of the result contains the new names.
- [x] A docs example: equations, renaming, LaTeX, and the alias-equation tip (`x ~ s.q`) for reading a renamed variable from the solution.

## Comments

- 2026-10-15: done. Keys may come from the uncompiled model (`model₊C₊q`, matched by stripping the root name) or the compiled one; an unknown key raises `ArgumentError`. Symbolics' `get_variables` returns derivative terms, which are mapped to their variable. Latexify and Symbolics are direct dependencies. Docs example: `examples/equations_latex.jl`, which also notes that a solution can't be indexed by a product of variables.
