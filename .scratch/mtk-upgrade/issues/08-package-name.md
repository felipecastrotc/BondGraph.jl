# Choose the package name

Type: grilling
Status: resolved
Blocked by: 01

## Question

What is the package called for v0.3 and registration: keep `BondGraph`, switch to `BondGraphToolkit` (as the README and docs already say), or another name, given BondGraphs.jl and the General registry's naming rules?

## Answer (grilling, 2026-10-15)

- **Name: `BondGraphToolkit`.** `BondGraph` fails all three AutoMerge similarity checks against the registered `BondGraphs`, and it collides with that package's exported `BondGraph` type (research/01). `BondGraphToolkit` passes the checks, follows SciML's `*Toolkit` precedent for "domain → MTK" packages, and is already the name in the README and docs. Soft risk: Python's BondGraphTools may draw a review comment (guideline §8); it is not an automatic failure.
- **Scope:** everything is renamed together and the UUID stays (`41e83443-…`, never registered):
  - module, `Project.toml` name, `src/BondGraphToolkit.jl`;
  - docs and tests;
  - GitHub repo `BondGraph.jl` → `BondGraphToolkit.jl` (AutoMerge requires the URL to end in `/<Name>.jl.git`).

  No `BondGraph` alias, since there are no users.
- **Timing:** the code rename lands on `feature/mtk-upgrade` with v0.3. The GitHub repo is renamed at the v0.3.0 release; ask the maintainer first, because it changes the docs URL.
