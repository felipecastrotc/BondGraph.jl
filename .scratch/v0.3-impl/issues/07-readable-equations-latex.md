# 07: Readable equations: `bg_equations` and `latexify`

**What to build:** A modeller gets the model's equations in terms of its states, renames variables for display (including inside derivatives, e.g. `D(s.q)` → `ẋ`, `m.p` → `p`), and exports them to LaTeX for a report or paper.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** ready-for-agent

- [ ] `bg_equations(sys; rename = Dict())` accepts a compiled or uncompiled model and returns MTK's full equations in terms of the states. The renames are applied everywhere, including inside `Differential`. The model is not modified.
- [ ] `latexify` is re-exported from Latexify, which is a direct dependency with a bound.
- [ ] Checks: renamed variables appear and the original names don't, including in derivatives; `latexify` of the result contains the new names.
- [ ] A docs example: equations, renaming, LaTeX, and the alias-equation tip (`x ~ s.q`) for reading a renamed variable from the solution.
