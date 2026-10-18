# Contributing to BondGraphToolkit.jl

Thank you for your interest. Bug reports, questions, ideas and pull requests are all welcome.

## Reporting a bug or asking a question

Open an [issue](https://github.com/felipecastrotc/BondGraph.jl/issues/new/choose). For a bug, include a minimal model that shows it, what you expected, what happened, and the output of `using Pkg; Pkg.status()` and `versioninfo()`.

A sign or equation that looks wrong is worth reporting too. The package follows the textbook conventions (Karnopp, Margolis & Rosenberg; Borutzky), and the [Concepts and Conventions](https://www.felipecastrotc.com/BondGraph.jl/stable/man/user_guide/) page states them.

## Proposing a change

1. Branch from `develop` as `feature/<kebab-name>`. `main` holds released code; `develop` is the integration branch.
2. Add a test for the behaviour you change, in `test/`, through the public API. The golden data in `test/regression/` comes from v0.2 and is never regenerated.
3. Run the tests and, if you touched the docs or examples, build the docs:
   ```sh
   julia --project -e 'using Pkg; Pkg.test()'
   julia --project=docs docs/make.jl
   ```
4. Open a pull request against `develop`.

Commit messages: an imperative, capitalized summary of at most 72 characters, without a trailing period (`Fix sign convention in Junction1`), and an optional body after a blank line.

Code style: the [Julia Base style guide](https://docs.julialang.org/en/v1/manual/style-guide/), with 4-space indentation, lowercase function names with words squashed where readable, `CamelCase` types and constructors, and `!` on mutating functions.

The [contributing page](https://www.felipecastrotc.com/BondGraph.jl/stable/contributing/) of the docs has the details on examples and figures.

By participating you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
