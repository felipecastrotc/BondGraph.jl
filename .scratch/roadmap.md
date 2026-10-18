# Roadmap: revive BondGraphToolkit.jl

Status: in progress (started 2026-10-14)

Steps 1 (commit pending work) and 2 (cleanup: deps, `[compat]`, CI, dead files) are done; see `git log` around `75c03c5`.

## Baseline (as of 2026-10-14)

- Last known-good stack: Julia 1.9.3, ModelingToolkit 8.71.2, Symbolics 5.8, SymbolicUtils 1.4 (local, gitignored `Manifest.toml`).
- With that Manifest, all 221 tests pass on Julia 1.10.
- Fresh resolve (no Manifest, what CI does) **fails**: a newer TimerOutputs breaks the old Groebner that Symbolics 5/MTK 8 depend on. CI is red until step 3 lands.
- `src/mtk.jl` overrides MTK methods. Julia 1.10+ refuses method overwriting during precompilation, so BondGraph can't be precompiled. This is the core of step 3.

## Step 3: upgrade to current ModelingToolkit

Branch: `feature/mtk-upgrade`. Destination: **BondGraph installs and passes tests on current MTK, Julia 1.10+, with fresh dependency resolution (CI green).**

- [x] `/research` (done: `.scratch/mtk-upgrade/research-mtk-internals.md`): what changed in MTK/Symbolics/SymbolicUtils from 8.71 to current for every internal `src/` imports (`src/mtk.jl`, `src/utils/*.jl`, `src/types.jl`). Output: a cited markdown note in the repo.
- [x] Regression safety net: golden-data tests for every example (`test/regression/`, report `.scratch/golden-tests/report.md`). Watch: N-DoF cancelling terms (0-junction merges non-adjacent elements) and Duffing's redundant state.
- [x] `/wayfinder`: map charted in `.scratch/mtk-upgrade/map.md` (10 decision tickets in `issues/`; the map is the source of truth from here). Original expected questions:
  - [x] Rebuild connection expansion on MTK's public API, no overrides (ADR 0001, tickets 03–05)
  - [x] Which MTK version to target, and the minimum Julia version: MTK `11.45` caret, Julia 1.10, CI 1.10 + 1 (ticket 09)
  - [x] Graph plotting: rebuilt on Graphviz with TikZ export; public API settled (ticket 10)
  - [x] Junction sign convention: v0.2 convention kept, confirmed against Karnopp/Borutzky (ticket 06)
  - [x] Package name: `BondGraphToolkit` (ticket 08)
- [ ] `/grill-with-docs` on the key decisions: produces `GLOSSARY.md` (effort, flow, port, junction, bond, ...) and ADRs in `docs/adr/`.
- [ ] TODO (future): SCAP causality as an optional diagnostic (ticket 10).
- [ ] `/to-spec` (done: `.scratch/mtk-upgrade/spec.md`) → `/to-tickets` (done: `.scratch/v0.3-impl/issues/01-10`) → `/implement` per ticket (done 2026-10-15, all 10 on `feature/mtk-upgrade`) (uses `tdd` and `code-review`).
- [x] Oriented ports (ticket 12's `p1` check): elements get `power`/`pout`, junction members `x`/`x => :pout`/`:a`/`:a => :in`; [ADR 0002](../docs/adr/0002-oriented-ports.md), `.scratch/oriented-ports/` (2026-10-18).
- [ ] `/diagnosing-bugs` for equation/sign mismatches against the literature (MSD N-DoF, DC motor).
- [ ] Fold the remaining `TODO.md` items (mTF, rename-by-Dict) into tickets, then delete `TODO.md`.
- [ ] TODO (deferred by maintainer): decide license and MTK vs ModelingToolkitBase dependency before registration/paper. AGPL is acceptable if needed; keep the current license until then. See `.scratch/mtk-upgrade/issues/07-dependency-license-stance.md`.
- [ ] Release: `release/0.3.0` → `main`, tag `v0.3.0`; rename the GitHub repo to `BondGraphToolkit.jl` (ask first; ticket 08).

## Step 4: documentation

Mostly after step 3: the examples must run on the new API.

Can start now:
- [x] Fix the README naming (depends on the package-name decision) (done 2026-10-18: README rewritten for v0.3)
- [x] Rewrite `docs/src/contributing.md` from `AGENTS.md` (Codex)
- [x] Docs/API coverage audit: `.scratch/codex/docs-audit/report.md` (Codex). Use it as the checklist for the items below.

After step 3:
- [ ] `/to-spec` for the docs scope, then `/to-tickets` → `/implement`:
  - [ ] User Guide (empty in `docs/make.jl`): concepts, elements, connections, simulation
  - [ ] API reference: docstrings for every export, `lib/public.md` and `lib/internals.md` complete; resolve exports with no local definition (`renamevars`, `addsubsys`, `latexify`), document `@oneport`
  - [ ] Fix `docs/make.jl` issues listed in the audit (`/tree/master/` link, cwd-relative paths, unguarded `rm`)
  - [ ] Examples: all run on the new API; finish `mechanical_msd_ndof.jl` (textbook `0 → sd`); decide on `tmp_dev_acoustic_carol.jl`
  - [ ] Fix the DC motor and 2-DoF TikZ drawings to match the textbook (ticket 06)
  - [ ] Changelog page (commented out in `docs/make.jl`)
  - [ ] Citing page with the thesis/paper reference (2026-10-18: cites the software from `CITATION.cff`; add the JOSS paper when published, see Publication)
- [x] Docs workflow green, deploy to `gh-pages` (v0.3.0 tag and `main` deployed, 2026-10-18)
- [ ] Optional: register in the Julia General registry

## Value gaps to reassess (2026-10-17)

From the assessment of whether the library adds value (user: record them and reassess later). Today it is mostly MTK components with textbook bond-graph conventions, readable equations and drawings. Its distinct strengths: elements are native MTK systems, and the MTKSL adapters. What would justify the bond-graph method over plain acausal MTK/MTKSL modelling is not built yet:

- [ ] Causality assignment (SCAP) as a diagnostic: `.scratch/v0.3-impl/issues/23-causality-diagnostic.md`
- [ ] Power and energy diagnostics: `.scratch/v0.3-impl/issues/24-power-energy-diagnostics.md`
- [ ] Multiport fields: `.scratch/v0.3-impl/issues/25-multiport-fields.md`
- [ ] Signals and modulated elements in drawings: `.scratch/v0.3-impl/issues/20-draw-signals.md`
- [x] Derivative-causality pitfalls (silent frozen pendulum): `.scratch/v0.3-impl/issues/21-derivative-causality-pitfalls.md` (documented; the example steers MTK with `state_priority`)

Positioning for the paper: engineering and multiphysics, textbook-faithful bond graphs as native MTK components interoperating with MTKSL, complementary to BondGraphs.jl (systems biology, graph object compiled to MTK, registered but effectively MTK 9 only; `.scratch/mtk-upgrade/research/01-bondgraphs-jl-comparison.md`). Other tools: BondGraphTools (Python), 20-sim (commercial).

## Publication: JOSS and Zenodo (after 0.3, decided 2026-10-18: wait for JOSS)

v0.3.0 is released (GitHub Release, `CITATION.cff`, community files, branch ruleset). Nothing needs a DOI before the JOSS paper, so Zenodo waits too.

Before submitting to JOSS:
- [ ] Decide the license and the MTK vs ModelingToolkitBase dependency (the deferred TODO above). JOSS accepts Apache-2.0 and AGPL, but changing it after publication is awkward.
- [ ] Register `BondGraphToolkit` in the Julia General registry (not required by JOSS, but reviewers will `pkg> add` it). Then add TagBot and CompatHelper workflows.
- [x] Add the ORCID to `CITATION.cff` (done 2026-10-18)
- [ ] Check JOSS's current requirements on its site: scope and "substantial scholarly effort", and its policy on AI-assisted code (likely a disclosure, and the history no longer has `Co-Authored-By` lines, so state it in the paper).
- [ ] Write `paper.md` and `paper.bib`: summary, statement of need, state of the field (BondGraphs.jl, BondGraphTools, 20-sim, plain MTK/MTKSL; positioning above), a short example and the references (Karnopp, Margolis & Rosenberg 2012; Borutzky 2011).
- [ ] Optional, before review: reassess the value gaps above (causality diagnostic, power/energy diagnostics) and the postponed example review (tickets 13–19). Reviewers will ask what the package adds over MTK.
- [ ] Submit at joss.theoj.org.

During review and on acceptance:
- [ ] Address the review in a release (e.g. v0.3.x or v0.4.0).
- [ ] Turn on the Zenodo GitHub integration (zenodo.org, sign in with GitHub, enable the repository) before that release, so it gets archived with a DOI. Later releases get one automatically; the concept DOI always points to the latest.
- [ ] Give JOSS the archive DOI and the version.
- [ ] After publication: add `preferred-citation` (the JOSS paper) to `CITATION.cff`, a JOSS and a DOI badge to the README, and the paper to `docs/src/citing.md`.
