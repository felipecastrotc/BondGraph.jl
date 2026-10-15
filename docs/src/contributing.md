# [Contributing](@id contribute)

Contributions and bug reports are welcome. Report bugs and propose changes in the [BondGraph.jl issue tracker](https://github.com/felipecastrotc/BondGraph.jl/issues).

## Working on a change

Create a `feature/<kebab-name>` branch from `develop`. The `develop` branch is for integration; `main` holds released code. Release branches use `release/<X.Y.Z>` from `develop`, and urgent fixes use `hotfix/<name>` from `main`.

Write commit summaries in the imperative, capitalized, at most 72 characters, with no trailing period (for example, `Fix sign convention in Junction1`). Put any explanatory body after a blank line.

Follow the [Julia Base style guide](https://docs.julialang.org/en/v1/manual/style-guide/): indent with four spaces, use `snake_case` for functions, and end mutating function names with `!`.

Run the tests from the repository root:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

## Documentation

Build the documentation from the repository root with `julia --project=docs docs/make.jl`.

Add examples as `examples/*.jl` Literate.jl sources. Start each file with `# # [Title](@id anchor)`; the title becomes its navigation entry. The example titled `Quick Example` is inserted into the getting-started page at `@quick-example`. Files prefixed `tmp_` are skipped.

For figures, keep TikZ sources in `docs/src/misc/tikz/*.tex`. Compile them to PDF, convert them to SVG and `-dark.svg` with `docs/src/misc/getsvg.sh`, and place the published SVGs in `docs/src/assets/tikz/`.

```@contents
Pages = ["contributing.md"]
Depth = 2
```
