# 10: Docs build green and API pages complete

**What to build:** A reader browses a docs site that builds cleanly and documents every v0.3 export with a docstring. Getting started uses the v0.3 API, and every element, port, macro and adapter has a runnable example.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 02, 03, 04, 05, 06, 07, 08, 09.

**Status:** done (feature/mtk-upgrade)

- [x] The public API page lists every export with its docstring; references to removed names are gone.
- [x] Getting started and the quick example use the v0.3 API.
- [x] A causality note: MTK chooses the states and reduces the index; a storage element in derivative causality gives a dummy-derivative variable; a causal conflict raises "structurally singular". Show both on small models.
- [x] Every element, port, macro and adapter has a runnable docs example (checklist against the export list).
- [x] The docs build passes locally and in CI.

## Comments

- 2026-10-15: done. Docs build locally (`julia --project=docs docs/make.jl`, exit 0; only warning: large plots fall back to SVG). `docs/Project.toml` rewritten for v0.3 (Documenter, Literate, Plots, OrdinaryDiffEq solvers, MTKSL, Graphviz_jll); `make.jl` fixes from the docs audit (paths from `@__DIR__`, guarded `rm`, `/tree/main/`, `edit_link`/`devbranch = "main"`). New page `man/user_guide.md` (Concepts and conventions) with the causality note run as `@example` blocks. Public API page lists every export; internals via `@autodocs Public = false`. Every export appears in runnable example code (checked against the export list); added a hand-written `PowerPort` element, translational/hydraulic adapter snippets and `bondgraph_dot` for that.
- Review of tickets 02-09 (subagent) fixed here: DC motor through `ElectricalAdapter` test added (ticket 06 criterion), `zero` local renamed, docs state that a plain `connect` to an element is a bond into it. Kept as is: the closed-loop signal test uses `De` and a closed form (stronger than comparing with the expression form).
- The old v0.2 docs Manifest was moved to `../Manifest-bondgraph-v0.2-docs.toml` (it blocked `Pkg.develop`).
