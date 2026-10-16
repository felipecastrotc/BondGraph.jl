# BondGraphToolkit.jl

Julia library for modeling and simulating dynamic systems with bond graphs, built on ModelingToolkit.jl (MTK 11). Module name `BondGraphToolkit` (v0.3 on branch `feature/mtk-upgrade`; v0.2 was `BondGraph` on MTK 8.71). Public MTK API only, no overrides of MTK internals (ADR 0001).

## Layout

- `src/BondGraphToolkit.jl`: module entry, exports, `PowerPort`, `@oneport`
- `src/elements.jl`: one-port elements (`Se`, `Sf`, `Resistance`, `Compliance`, `Inertia`) and the mechanical shorthands `Mass`, `Spring`, `Damper`
- `src/junctions.jl`: the signed 0-/1-junction component (`Junction0`, `Junction1`)
- `src/two_ports.jl`: `TF`, `GY` (built with `@twoport`, defined in the module file)
- `src/adapters.jl`: the five MTK standard-library adapters
- `src/graph.jl`: `bondgraph` (the model's bond graph as data, shown as text), `bgplot` (SVG/TikZ/DOT drawing, laid out with Graphviz_jll), `savebg` (export by file extension)
- `src/equations.jl`: `bg_equations` (readable, renameable equations; `latexify` is re-exported)
- `test/runtests.jl` + `test/test_*.jl`: behaviour tests through the public API; `test/regression/`: golden data from v0.2 (never regenerate)
- `examples/*.jl`: Literate.jl sources, rendered into the docs
- `docs/`: Documenter.jl site, deployed to `gh-pages`; has its own `docs/Project.toml`

## Commands

```sh
julia --project -e 'using Pkg; Pkg.test()'      # tests (CI: Julia 1.10 and 1, no Manifest)
julia --project=docs docs/make.jl               # build docs (also tries to deploy; harmless locally)
```

## Conventions

- Style: Julia Base style guide, 4-space indent, `snake_case` functions, `!` on mutating functions.
- Examples: first line is `# # [Title](@id anchor)`; the title becomes the docs nav entry. The example titled `Quick Example` is spliced into `docs/src/getting_started.md` at `@quick-example`. Files prefixed `tmp_` are skipped.
- Figures: TikZ sources in `docs/src/misc/tikz/*.tex` → PDF → SVG (+ `-dark.svg`) via `docs/src/misc/getsvg.sh`; published copies go in `docs/src/assets/tikz/`.
- `tmp_dev_*` files are personal scratch, gitignored. Never commit them.
- Dependencies: runtime deps in `Project.toml` `[deps]` with a `[compat]` bound each; test-only deps go in `[extras]`/`[targets]`; example/plotting deps go in `docs/Project.toml`. `Manifest.toml` is gitignored.

## Git workflow

git-flow naming (prefixes recorded in `.git/config`; the `git flow` CLI isn't required, plain git works):

- `main` = released code, `develop` = integration branch.
- New work: `feature/<kebab-name>` off `develop`, merged back with `--no-ff`.
- Releases: `release/<X.Y.Z>` off `develop`; bump `version` in `Project.toml`, merge to `main` and `develop`, tag `vX.Y.Z` (annotated, short description).
- Urgent fixes: `hotfix/<name>` off `main`, merged to both.
- Versions follow SemVer (required by the Julia General registry). Legacy tags `vMTK-v8.71.2*` are dependency-pinned hotfixes; don't repeat that scheme.
- Commit messages: imperative summary ≤72 chars, capitalized, no trailing period ("Fix sign convention in Junction1"), optional body after a blank line.
- Don't push, tag, or merge to `main` without asking.

## Agent skills

### Issue tracker

Issues live as local markdown files under `.scratch/<feature>/`. See `docs/agents/issue-tracker.md`.

### Triage labels

Default five-role vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`), recorded as `Status:` lines. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `GLOSSARY.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.

Bond-graph theory (sign conventions, element laws, causality): consult `docs/literature/README.md` and the cited book pages before answering or changing conventions. Keep what a book states separate from inference and design choices.
