# 01: Tracer bullet: BondGraphToolkit on MTK 11 with the RC example

**What to build:** A modeller installs `BondGraphToolkit` with fresh dependency resolution on Julia 1.10+ and MTK 11. They build the RC circuit from `Se`/`Sf`, `Resistance`, `Compliance` and a signed junction, then simulate it with `mtkcompile` → `ODEProblem` → `solve`. The results match v0.2's golden data. This is the core every other ticket builds on.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`. Prior art: the prototype on branch `prototype/connection-mechanism` (`BGProto.jl` and its check script) and ADR 0001.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Branch `feature/mtk-upgrade` off `develop`.
- [ ] Package, module, main source file, README and tests renamed to `BondGraphToolkit`; UUID unchanged; no alias.
- [ ] `[compat]`: `ModelingToolkit = "11.45"`, `julia = "1.10"`. Direct dependencies only for what the code imports, each bounded. `@reexport using ModelingToolkit` is kept.
- [ ] Deleted: the MTK overrides, the bond metadata types, the connection and graph utilities, `simplifysys`/`reducedobs`/`substitute_dict`, the `structural_simplify` export, `renamevars`/`addsubsys`/`latexify` dead exports, `Dq`, `isindependent`, `It`/DomainSets, `generate_graph` with Plots/GraphRecipes, `src/mbg.jl`, `mGY`/`mTF`, the old `Mass`/`Spring`/`Damper`. Their unit tests are deleted too.
- [ ] `PowerPort` (`e`, and `f` marked `Flow`); MTK's `t_nounits`/`D_nounits` are used, with no package-owned `t`/`D`.
- [ ] One signed junction component behind `Junction0`/`Junction1`. Members: `x` (+1, into the junction), `[-1, x]` (out), and free ports `:name`/`[-1, :name]`. Equations as in the spec: shared variable unsigned, signed sum.
- [ ] `@oneport` adds `power` and wraps the laws in a `System`. `Se(E)`, `Sf(F)`, `Resistance(R)`, `Compliance(C, q)` and `Inertia(I, p)` are written with it, with laws in bond variables and docstrings showing the conventional drawing.
- [ ] The RC example is rewritten for v0.3 and reproduces its golden CSV. State names are remapped only in the regression helpers.
- [ ] Physics checks through the public API:
  - [ ] Junction power balance `Σ bᵢ eᵢ fᵢ = 0` along a trajectory, at a 0- and at a 1-junction.
  - [ ] A Resistance drawn reversed delivers power (`q = e^{+t}` in the RC case).
- [ ] The regression tests of not-yet-ported examples are skipped by name, with a pointer to their ticket.
- [ ] CI is green on Julia `1.10` and `1` with no Manifest.
