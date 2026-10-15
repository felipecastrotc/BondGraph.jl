# Roadmap: revive BondGraph.jl

Status: in progress (started 2026-10-14)

Steps 1 (commit pending work) and 2 (cleanup: deps, `[compat]`, CI, dead files) are done; see `git log` around `75c03c5`.

## Baseline (as of 2026-10-14)

- Last known-good stack: Julia 1.9.3, ModelingToolkit 8.71.2, Symbolics 5.8, SymbolicUtils 1.4 (local, gitignored `Manifest.toml`).
- With that Manifest, all 221 tests pass on Julia 1.10.
- Fresh resolve (no Manifest, what CI does) **fails**: a newer TimerOutputs breaks the old Groebner that Symbolics 5/MTK 8 depend on. CI is red until step 3 lands.
- `src/mtk.jl` overrides MTK methods. Julia 1.10+ refuses method overwriting during precompilation, so BondGraph can't be precompiled. This is the core of step 3.

## Step 3: upgrade to current ModelingToolkit

Branch: `feature/mtk-upgrade`. Destination: **BondGraph installs and passes tests on current MTK, Julia 1.10+, with fresh dependency resolution (CI green).**

- [ ] `/research`: what changed in MTK/Symbolics/SymbolicUtils from 8.71 to current for every internal `src/` imports (`src/mtk.jl`, `src/utils/*.jl`, `src/types.jl`). Output: a cited markdown note in the repo.
- [ ] `/wayfinder`: chart the upgrade as a map in `.scratch/mtk-upgrade/map.md` with decision tickets. Expected questions:
  - [ ] Keep overriding MTK internals, or rebuild connection expansion on MTK's public API (the big one)
  - [ ] Which MTK version to target, and the minimum Julia version
  - [ ] Turn `generate_graph` (Plots/GraphRecipes) into a package extension so the core stays light
  - [ ] Junction sign convention (see `TODO.md`: DC motor matches literature only with the sign disabled)
  - [ ] Package name: `BondGraph` vs `BondGraphToolkit` (`BondGraphs.jl` already exists in the General registry)
- [ ] `/grill-with-docs` on the key decisions: produces `GLOSSARY.md` (effort, flow, port, junction, bond, ...) and ADRs in `docs/adr/`.
- [ ] `/to-spec` → `/to-tickets` → `/implement` per ticket (uses `tdd` and `code-review`).
- [ ] `/diagnosing-bugs` for equation/sign mismatches against the literature (MSD N-DoF, DC motor).
- [ ] Fold the remaining `TODO.md` items (mTF, rename-by-Dict) into tickets, then delete `TODO.md`.
- [ ] Release: `release/0.3.0` → `main`, tag `v0.3.0`.

## Step 4: documentation

Mostly after step 3: the examples must run on the new API.

Can start now:
- [ ] Fix the README naming (depends on the package-name decision)
- [x] Rewrite `docs/src/contributing.md` from `AGENTS.md` (Codex)
- [x] Docs/API coverage audit: `.scratch/codex/docs-audit/report.md` (Codex). Use it as the checklist for the items below.

After step 3:
- [ ] `/to-spec` for the docs scope, then `/to-tickets` → `/implement`:
  - [ ] User Guide (empty in `docs/make.jl`): concepts, elements, connections, simulation
  - [ ] API reference: docstrings for every export, `lib/public.md` and `lib/internals.md` complete; resolve exports with no local definition (`renamevars`, `addsubsys`, `latexify`), document `@oneport`
  - [ ] Fix `docs/make.jl` issues listed in the audit (`/tree/master/` link, cwd-relative paths, unguarded `rm`)
  - [ ] Examples: all run on the new API; finish `mechanical_msd_ndof.jl`; decide on `tmp_dev_acoustic_carol.jl`
  - [ ] Changelog page (commented out in `docs/make.jl`)
  - [ ] Citing page with the thesis/paper reference
- [ ] Docs workflow green, deploy to `gh-pages`
- [ ] Optional: register in the Julia General registry
