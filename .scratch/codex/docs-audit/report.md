# Documentation and API coverage audit

Scope: explicit `export` statements in `src/BondGraph.jl`, local source docstrings, and the checked-in documentation. `@reexport using ModelingToolkit` may expose additional upstream symbols; those are outside this explicit-export inventory. “Yes” means a local docstring is attached to the definition, not that the example or API was executed. No docstring here is merely an empty placeholder. The docs build was not run, as requested.

## Explicit exports

| Symbol | Docstring? | In `lib/public.md`? | What the source does |
| --- | --- | --- | --- |
| `Power` | yes | yes | Creates an effort/flow connector with bond-graph metadata. |
| `Se` | yes | yes | Creates a source that fixes a port's effort to an expression. |
| `Sf` | yes | yes | Creates a source that fixes a port's flow to an expression. |
| `Dq` | no | no | Creates an ODE system with `D(q) ~ f`; its `x` keyword is unused. |
| `Junction1` | yes | yes | Connects signed subsystems through a one-junction port. |
| `Junction0` | yes | yes | Connects signed subsystems through a zero-junction port. |
| `mGY` | yes | yes | Builds a two-port gyrator and its connection equations. |
| `mTF` | yes | yes | Builds a two-port transformer and its connection equations. |
| `Mass` | yes | yes | Builds an inertial element with `D(flow) ~ effort / I`. |
| `Spring` | yes | yes | Builds a compliant element with displacement as a state. |
| `Damper` | yes | yes | Builds a resistive element with `effort ~ flow * R`. |
| `isindependent` | yes | yes | Tests whether a symbolic value is independent of the model's time variable. |
| `generate_graph` | yes | yes | Plots a graph of effort or flow connections. |
| `simplifysys` | yes | yes | Expands connections, structurally simplifies, then substitutes observed variables. |
| `renamevars` | no | no | Has no active local definition; a proposed renaming implementation is commented out in `src/mtk.jl`. |
| `latexify` | no local docstring | no | Has no active local definition; examples use it to format equations as LaTeX. Its binding appears intended to come from a dependency. |
| `structural_simplify` | no local docstring | no | Imported from ModelingToolkit; simplifies the structure of a symbolic system. |
| `addsubsys` | no | no | Has no active local definition or use in this repository. |
| `@oneport` | no | no | Rewrites a one-port function body to create `Power`, collect equations, and compose an `ODESystem`. |

The public page also lists `BondGraph` (the module itself) and `BondGraph.GenericDamper` (documented but not exported). `reducedobs` has a docstring and appears in the module docstring's prose export list, but is absent from the actual `export` statements and public page. `Se` and `Sf` have substantive docstrings with “TODO: Improve” notes; they are counted as yes, not stubs. For `renamevars`, `latexify`, and `addsubsys`, actual runtime binding was not verified because loading the dependency stack is outside this audit.

## Internals

`lib/internals.md` uses `@autodocs` for four source pages. The following non-exported functions have docstrings there and are selected by those directives:

| Source page | Documented non-exports |
| --- | --- |
| `utils/bg_utils.jl` | `equalityeqs`, `sumvar`, `flatinput` |
| `utils/connection_utils.jl` | `check_bg_con`, `get_bg_connection_set!`, `add_idx`, `gen_tp_con!` |
| `utils/graph.jl` | `generate_bg_eqs!`, `csets2dict`, `csets2adjmtx`, `adjmtx2eqs`, `retrieve_node_details`, `gen_con_name` |
| `utils/system_utils.jl` | `reducedobs`, `substitute_dict`, `substitute_observed` |

Important missing internals: `get_var` in `utils/graph.jl` has no docstring, so the graph-page `@autodocs` omits it; `set_bg_metadata`, `get_bg_junction`, and `update_mtk_con` in `src/types.jl` define the bond-graph metadata layer but have neither docstrings nor an internals-page entry; the active ModelingToolkit overrides `expand_connections`, `generate_connection_set`, `split_assign`, and `single_named_expr` in `src/mtk.jl` are also undocumented and not included. `@oneport` is exported but absent from both API pages. `GenericDamper` is documented on the public page despite being non-exported. The `@autodocs` source filter has not been checked by a Documenter build.

## Empty and placeholder pages

| Page or navigation entry | Current contents |
| --- | --- |
| User Guide (`docs/make.jl:132-133`) | Empty navigation group. `docs/src/man/user_guide.md` exists but is unlinked: a short introduction followed by unfinished headings and an incomplete “The” paragraph. |
| Changelog (`docs/make.jl:142`) | Commented-out navigation entry; no `docs/src/changelog.md` exists. |
| Citing (`docs/src/citing.md`) | Heading, request to cite a paper, `TODO: add the paper`, and an article link; no complete citation. |
| Getting started (`docs/src/getting_started.md`) | Installation steps using `pkg> add BondGraphToolkit.jl`, `TODO: add instalation latest version`, and an `@quick-example` marker. `docs/make.jl` replaces that marker with the generated RC example in `getting_started_auto.md`. |

## `docs/make.jl` issues and suggested fixes

| Location | Issue | Suggested fix |
| --- | --- | --- |
| 40 | Example source link uses `/tree/master/`; the branch convention says `main`. | Change link to `/tree/main/examples/`. |
| 85-91, 114 | `get_nav` opens `./examples/` relative to the process working directory, while generation uses `@__DIR__`. | Resolve example paths from `@__DIR__` consistently. |
| 116-119 | `electrical_rc.jl` is titled “Quick Example”; after generating it, the script inserts it into getting started and unconditionally removes `docs/src/examples/electrical_rc.md`. A renamed or missing generated file would make this `rm` fail. | Derive the page path from the selected quick example and remove it only if it exists, or skip generating a standalone page for that example. |
| 59-81, 116-117 | `insert_quick_example(quick_example)` assumes exactly one matching title, a heading in its generated page, and an `@quick-example` marker. | Validate these inputs and report a clear error when the example or marker is missing. |
| 132-133, 142 | Empty User Guide group and commented Changelog leave authored or planned pages outside navigation. | Link a completed User Guide and add a real Changelog page when written; otherwise remove the empty group. |
| 25-29, 151-152 | Build script removes `docs/build` and `getting_started_auto.md` using hardcoded paths and unconditional cleanup. | Use `joinpath(@__DIR__, ...)` for build paths and guard cleanup of generated files. |
| 8-9, 155 | Local build calls `Pkg.develop`/`Pkg.instantiate` and reaches `deploydocs` in the same script. | Separate dependency setup and deploy into explicit CI steps once the docs workflow is repaired. |

## Naming occurrences in README and docs

The package name is undecided. The following are all occurrences of `BondGraph.jl` or `BondGraphToolkit` in `README.md` and checked-in docs, including URLs and repeated mentions on one line. No naming edits were made outside the requested contributing-page rewrite.

| File:line | Occurrence |
| --- | --- |
| `README.md:9` | `BondGraph.jl` in documentation URL |
| `README.md:16` | `BondGraphToolkit.jl` twice |
| `docs/make.jl:40` | `BondGraph.jl` in example URL |
| `docs/make.jl:126` | `BondGraph.jl` site name |
| `docs/make.jl:139` | `BondGraphToolkit` in Citing navigation title |
| `docs/make.jl:155` | `BondGraph.jl` in deploy repository URL |
| `docs/src/citing.md:1` | `BondGraphToolkit` in heading |
| `docs/src/citing.md:3` | `BondGraphToolkit.jl` |
| `docs/src/contributing.md:3` | `BondGraph.jl` in link text and issue URL |
| `docs/src/getting_started.md:14` | `BondGraphToolkit.jl` in installation command |
| `docs/src/index.md:5` | `BondGraphToolkit.jl` twice |
| `docs/src/index.md:33` | `BondGraphToolkit.jl` in installation command |
| `docs/src/lib/internals.md:3` | `BondGraph.jl` |
| `docs/src/lib/public.md:4` | `BondGraph.jl` |
| `docs/src/license.md:3` | `BondGraphToolkit` |
| `docs/src/license.md:22` | `BondGraph.jl` in license URL |
| `docs/src/license.md:24` | `BondGraphToolkit` |

STATUS: DONE
