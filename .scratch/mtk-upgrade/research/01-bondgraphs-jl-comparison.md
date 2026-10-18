# Research 01: BondGraphs.jl today, and what it means for our name

Ticket: [issues/01-bondgraphs-jl-comparison.md](../issues/01-bondgraphs-jl-comparison.md). Researched 2026-10-14.

## Summary

1. **BondGraphs.jl** (General name `BondGraphs`, UUID `d1a4f103-cd9c-465b-aef3-7f12ef11bc98`, MIT) is by Joshua Forrest (Univ. of Melbourne; with M. Pan, M. P. H. Stumpf, V. Rajagopal). Its paper is *Bioinformatics* 39(10) btad578 (2023). Latest release: v0.5.1 (2025-12-19). Active unreleased rewrite on `dev` (last commit 2026-04-08).
2. Its model is a **Graphs.jl graph object** (`BondGraph`, `Component(:R)`, `connect!`), not an MTK system. MTK is a backend: `ODESystem(bg)` / `simulate(bg, tspan)` compile it at the end. Catalyst is used to convert `ReactionSystem`s into bond graphs. Its focus is **systems biology / biochemistry**.
3. Custom elements are **Dict-schema library entries** (`addlibrary!`). Built-ins are R, C, I, Se, Sf and a constant-ratio TF, plus the biochemical Ce/Re/SCe. There is **no GY, no modulated TF/GY**, and controls are functions of `t` only. The authors call the custom-component API "unintuitive" (issue #68).
4. Overlap with us: both are bond-graph libraries that emit MTK equations. Difference: we make **every element and model an ordinary MTK system**, composed with `@named`/`compose`. They build a separate graph type and convert it to MTK. They target biochemistry; we target engineering and multiphysics with nonlinear R/C and mTF/mGY.
5. In practice it resolves only on **MTK 9**. Its compat claims MTK 9–11, but `Catalyst = "14, 15"` (both need MTK 9) and `SymbolicUtils = "3"` (MTK 11 needs v4) rule the others out.
6. **Name `BondGraph` fails all three AutoMerge distance checks against `BondGraphs`**: DL = 1, lowercase DL = 1, visual distance 2.27 ≤ 2.5. It is also the textbook case of naming rule 7 ("Websocket" vs "WebSockets"). On top of that, `BondGraphs` *exports a type called `BondGraph`*, so the two names would collide in user code. Treat `BondGraph` as not registrable.
7. `BondGraphToolkit`, `BondGraphsToolkit`, `BondGraphModelingToolkit`, `ModelingToolkitBondGraphs`, `BondGraphKit` and `BondGraphMTK` all **pass** the automated checks against all 12,675 non-JLL General names. `*Toolkit` has SciML precedent (`CellMLToolkit`, `SBMLToolkit`). `MTK` is an acronym and so falls under the rule 1 "avoid jargon" risk.
8. Soft risk for `BondGraphToolkit`: the Python package **BondGraphTools** ("A Toolkit for modelling multi-physics systems") is in the same domain, which may draw a rule-8 comment during the 3-day review. This is judgement, not an automated check.
9. Other bond-graph Julia packages: only `BondGraphs` is registered (plus the unrelated `CrystalBonds`). `SystemBondGraphs.jl` (TRACER-LULab, MTK-based) exists on GitHub but is unregistered and has been inactive since 2024-06.

## 1. BondGraphs.jl: identity and status

| Fact | Value | Source |
|---|---|---|
| Registered name / UUID | `BondGraphs`, `d1a4f103-cd9c-465b-aef3-7f12ef11bc98` | [General `B/BondGraphs/Package.toml`](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/B/BondGraphs/Package.toml); [Registry.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/Registry.toml) |
| Repo | https://github.com/jedforrest/BondGraphs.jl (default branch `master`, created 2021-08-18, 20 stars, 8 open issues as of 2026-10-14) | `gh api repos/jedforrest/BondGraphs.jl` |
| License | MIT, "Copyright (c) 2021 Jedforrest" | [LICENSE@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/LICENSE) |
| Registered versions | 0.3.0, 0.3.1, 0.4.0, 0.5.0, 0.5.1 | [Versions.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/B/BondGraphs/Versions.toml) |
| Release dates | v0.3.0 2023-03-05, v0.4.0 2023-03-24, v0.5.0 2024-12-19, **v0.5.1 2025-12-19** ("Refactoring source code and tests in preparation for a major refactor", PR #85) | [releases](https://github.com/jedforrest/BondGraphs.jl/releases); [commit eeac405](https://github.com/jedforrest/BondGraphs.jl/commit/eeac40590943cb71edfd2dee1a9c21f491baeb00) |
| Ongoing work | `dev` / `simulation-tests` branches, last commit 2026-04-08 ("WIP simulation tests"). New `BondGraphElement`, `Port`, MetaGraphsNext backend, function-style element constructors in `libraries/electrical.jl` | [dev@78c4e5c](https://github.com/jedforrest/BondGraphs.jl/tree/78c4e5c2e9f4ebdc8520a51e4e1431190fcd0452) |
| Open roadmap issues | #83 bond-graph ontology (CellML/OWL), #86 faster system generation, #87 MetaGraphsNext backend, #88 "New `struct` for defining component equations" (replace the Dict libraries), #90 SCAP causality assignment, #68 "Adding custom BG components" ("unintuitive and not very flexible") | [issues](https://github.com/jedforrest/BondGraphs.jl/issues) |
| Contributors | jedforrest (196 commits), mic-pan (55) | `gh api repos/jedforrest/BondGraphs.jl/contributors` |
| Docs | Documenter site with stable = v0.5.1 and versions v0.3–v0.5 | https://jedforrest.github.io/BondGraphs.jl/stable ; [gh-pages versions.js](https://github.com/jedforrest/BondGraphs.jl/blob/gh-pages/versions.js) |

**Paper.** J. Forrest, V. Rajagopal, M. P. H. Stumpf, M. Pan, "BondGraphs.jl: composable energy-based modelling in systems biology", *Bioinformatics* 39(10), btad578. It is an Applications Note, published 2023-09-19, CC BY 4.0, DOI [10.1093/bioinformatics/btad578](https://doi.org/10.1093/bioinformatics/btad578). All authors are at the University of Melbourne (School of BioSciences, School of Mathematics and Statistics, Biomedical Engineering); Pan and Stumpf are corresponding authors. The abstract describes models as "mutable, hierarchical, multiscale, and multiphysics" and compatible with the Julia modelling ecosystem. Source: [OUP article page](https://academic.oup.com/bioinformatics/article/39/10/btad578/7276395); the README "Cite" section also gives this reference ([README@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/README.md)).

**Dependencies and effective MTK version.**
- [Project.toml@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/Project.toml) depends on Catalyst, ModelingToolkit, DifferentialEquations, Graphs, GraphMakie, GraphRecipes, Latexify, SymbolicUtils, Setfield and OrderedCollections. Its compat is `ModelingToolkit = "9, 10, 11"`, `Catalyst = "14, 15"`, `SymbolicUtils = "3"`, `julia = "1"`. The registry records the same ([Compat.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/B/BondGraphs/Compat.toml)).
- Catalyst 14.x and 15.x both require `ModelingToolkit` 9.x ([Catalyst Compat.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/C/Catalyst/Compat.toml), e.g. `["15.0.9 - 15"] ModelingToolkit = "9.73.0 - 9"`). Catalyst 16 moved to `ModelingToolkitBase`.
- MTK 11.x requires `SymbolicUtils` 4.x ([MTK Compat.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/M/ModelingToolkit/Compat.toml)).
- **Inference:** BondGraphs 0.5.1 can only resolve with MTK 9. I did not run a Pkg resolve, so this rests on the registry compat entries above.
- The code uses `ODESystem`, `structural_simplify` and `@set! model.eqs`/`model.substitutions` ([systems.jl@v0.5.1 L175-191](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L175-L191)).

## 2. Scope and API (v0.5.1)

- **Model = a graph object, not an MTK system.** `BondGraph` is its own struct on Graphs.jl. You fill it imperatively: `add_node!`, `connect!`, `swap!`, `insert_node!`, `merge_nodes!`, `simplify_junctions!`, `expose`. You read it with `constitutive_relations`. You run it with `simulate`. All of these are exported in [src/BondGraphs.jl@v0.5.1 L41-L86](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/BondGraphs.jl#L41-L86). Graphs.jl functions such as `incidence_matrix(model)` work on it ([gettingstarted.md](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/docs/src/gettingstarted.md)).
- **How it becomes MTK.**
  - `ModelingToolkit.ODESystem(n::AbstractNode)` builds one `ODESystem` per node. Each port is an MTK `@connector MTKPort` with `E(t)` and `F(t) [connect = Flow]`.
  - `ODESystem(bg::BondGraph)` turns every bond into `connect(src_port, dst_port)`, composes the result and, by default, runs `structural_simplify`.
  - `simulate` wraps `ODEProblem` + `solve` ([systems.jl@v0.5.1 L83-L209](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L83-L209)).
  - The docs say: "Under the hood, our `simulate` function is converting our bond graph into an `ModelingToolkit.ODESystem`. We can chose instead to create an `ODESystem` directly" ([gettingstarted.md L86](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/docs/src/gettingstarted.md)).
- **Hierarchy.** A `BondGraphNode` wraps a sub-bond-graph. Its external ports come from `SourceSensor` nodes ([systems.jl L148-L165](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L148-L165); [AbstractNode.jl L61-L74](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/structures/AbstractNode.jl#L61-L74)).
- **Catalyst.** `BondGraph(rs::ReactionSystem; chemostats=[])` converts a reaction network into Ce/Re/TF/junction structure ([catalyst.jl@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/catalyst.jl)). This is the bridge in that direction only: Catalyst → bond graph.
- **Element library** ([standard.jl@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/libraries/standard.jl), [biochemical.jl](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/libraries/biochemical.jl)):
  - Elements are `Dict`s with `:description`, `:numports`, `:variables` (`:parameters`/`:globals`/`:states`/`:controls`) and `:equations`, written over the generic port arrays `E(t)[1:N]`, `F(t)[1:N]`.
  - Standard elements: linear `:R`, `:C`, `:I`, `:Se`, `:Sf`, and `:TF` with a constant parameter `n` (`e₂ = n*e₁`, `f₁ = -n*f₂`).
  - Biochemical elements: `:Ce`, `:ce`, `:Re`, `:re`, `:SCe`.
  - **No `GY`, no modulated TF/GY.** A grep for `GY|gyrator|modulat` over `src` and `docs` at v0.5.1 returns nothing.
  - `:controls` are "time-dependent parameters, can accept an arbitrary julia function", set as `Is.fs = t -> sin(2t)` ([standard.jl L14](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/libraries/standard.jl#L14); [gettingstarted.md L98-L103](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/docs/src/gettingstarted.md)).
  - UNVERIFIED: nothing documented lets one element's ratio depend on another element's state, i.e. true mTF/mGY. The schema only refers to the element's own port variables, its own variables and globals.
- **Custom elements.** You build a Dict and register it into the global `DEFAULT_LIBRARY` with `addlibrary!(Dict(:ReMM => ReMM))`. Then you instantiate it with `Component(:ReMM, "MM1"; r1=100)`. The docstring warns: "This library is likely to change in the future. Do not rely on this schema too much." Sources: [libraryfunctions.jl@v0.5.1](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/libraries/libraryfunctions.jl), [examples.md L330-L400](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/docs/src/examples.md). Issue [#68](https://github.com/jedforrest/BondGraphs.jl/issues/68) calls this "unintuitive and not very flexible".
- **Direction of the rewrite (unreleased).** `dev` adds function-style constructors such as `resistor(; R=1, name)` → `DissipatorElement([V ~ R * I], [V], [I]; name)` ([dev electrical.jl](https://github.com/jedforrest/BondGraphs.jl/blob/78c4e5c2e9f4ebdc8520a51e4e1431190fcd0452/src/libraries/electrical.jl)). This moves closer to MTK idioms, but the model is still a separate `BondGraph`/`BondGraphElement` object compiled through `ODESystem(...)` ([dev systems.jl](https://github.com/jedforrest/BondGraphs.jl/blob/78c4e5c2e9f4ebdc8520a51e4e1431190fcd0452/src/systems.jl)).

## 3. Overlap and difference with BondGraphToolkit.jl's goal

Our goal, from [map.md](../map.md): bond-graph models *are* ordinary MTK systems. They compose and couple with any MTK model, and custom nonlinear R/C and mGY/mTF stay easy to write.

| Axis | BondGraphs.jl v0.5.1 | BondGraphToolkit.jl (this repo, v0.2.1 @ `7638f46`) |
|---|---|---|
| Model object | Own `BondGraph` graph type (Graphs.jl), mutated with `add_node!`/`connect!`; MTK only at compile time | Each element is an MTK `ODESystem` (`compose(ODESystem(eqs,t,…), power)`, [src/bg.jl L406-L441](../../../src/bg.jl)); a junction composes MTK systems |
| Building style | Imperative graph editing (`swap!`, `insert_node!`, …), graph algorithms | MTK-style `@named x = Element(...)` + junction constructors |
| Custom element | Dict schema + global `addlibrary!` | Plain Julia function returning an MTK system, or `@oneport function Spring3(; name, k, x) … end` ([examples/mechanical_duffing.jl L55](../../../examples/mechanical_duffing.jl)) |
| Nonlinear R/C | Possible via Dict equations (e.g. biochemical `Ce` uses `log`) | First-class (`@oneport`, `GenericDamper(expr)`) |
| Modulated two-ports | None (constant `TF` only, no `GY`) | `mGY(...; g)` / `mTF(...; r)` accept symbolic expressions of system variables ([src/bg.jl L254-L330](../../../src/bg.jl)) |
| Interop with arbitrary MTK models | Indirect: convert to `ODESystem` first. Ports are named `p1…pN` with fields `E`/`F` (UNVERIFIED: no documented example of connecting to non-BG MTK components) | The model already *is* an MTK system (the coupling mechanism is an open question in map.md) |
| Domain focus | Systems biology / biochemistry, Catalyst import, thermodynamics (paper title) | Engineering / multiphysics (mechanical, electrical, DC motor, Duffing examples) |
| MTK version | MTK 9 in practice (see §1) | `~8.71` today; v0.3 targets MTK 11 |
| Causality | Not assigned; SCAP planned (issue #90) | Left to MTK index reduction (map.md "Not yet specified") |
| License | MIT | Apache-2.0 ([LICENSE](../../../LICENSE)) |

### Side by side: RC circuit

BondGraphs.jl. Source: [docs/src/examples.md@v0.5.1, "Simple Electric Circuit"](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/docs/src/examples.md); the current source is omitted here.

```julia
using BondGraphs
model = BondGraph("RC Circuit")
C = Component(:C)
R = Component(:R)
kvl = EqualEffort()                     # 0-junction
add_node!(model, [C, R, kvl])
connect!(model, R, kvl)
connect!(model, C, kvl)
C.C = 1
R.R = 2
sol = simulate(model, (0., 10.); u0 = [1])   # u0 for C.q
```

BondGraphToolkit.jl. Source: [examples/electrical_rc.jl L27-L48 @ 7638f46](../../../examples/electrical_rc.jl).

```julia
using BondGraph, DifferentialEquations
@named C = Spring(k = 1 / 1.0)       # compliance
@named R = Damper(c = 1.0)           # resistance
@named model = Junction0(R, C)       # an MTK system
@named sys = simplifysys(model)
prob = ODEProblem(sys, [R.power.f => 0.0, C.q => 1.0], (0, 10))
sol = solve(prob)
```

What the comparison shows:
- In BondGraphs.jl you build a graph and call `simulate`. In BondGraphToolkit.jl the junction call returns an MTK system, and you use the standard `ODEProblem`/`solve`.
- BondGraphs.jl indexes initial conditions positionally (`u0 = [1]`). BondGraphToolkit.jl uses MTK symbolic maps.

## 4. Naming: General registry rules and AutoMerge checks

**Rules (General [NAMING_GUIDELINES.md](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/NAMING_GUIDELINES.md)):**
- §1: avoid jargon and acronyms.
- §3: pluralize when the package mainly provides a new type, e.g. `DataFrames` provides `DataFrame`.
- §4: prefer clarity, with names of at least 5 characters.
- §7: "Avoid naming a package closely to an existing package. `Websocket` is too close to `WebSockets`… Rather use a new name… that is descriptive of what distinguishes the package."
- §8: avoid a distinctive name already used by a well-known unrelated project.
- §9: UpperCamelCase, ASCII only.

The General [README](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/README.md) adds these points:
- New packages wait 3 days so the community can check, among other things, that the name is not typosquatting and fits the naming guidelines.
- Earlier exceptions do not count as precedent.
- Renaming after registration means re-registering under a new UUID.

**AutoMerge name checks** ([RegistryCI AutoMerge `guidelines.jl` @ 0cb8b76](https://github.com/JuliaRegistries/RegistryCI.jl/blob/0cb8b7615cbd5db46b9652a1d7fc048b0aea25d6/AutoMerge/src/guidelines.jl), tag AutoMerge-v1.1.1):
- Name length ≥ 5 (L243).
- Valid identifier and ASCII.
- `^[A-Z][A-Za-z0-9]*?$` capitalization (L703).
- No `julia`, `Ju` or `jl` in the name.
- An exact or case-insensitive match is a hard fail (`meets_name_match_check`).
- Repo URL must end in `/<Name>.jl.git` (L1317). A rename therefore also means renaming the GitHub repo.
- **Distance check** `meets_distance_check` (L565-L632) compares the name against every non-JLL name. It fails if any of these holds:
  - Damerau–Levenshtein ≤ 2 (`StringDistances.DamerauLevenshtein`, L292).
  - Lowercased DL ≤ 1.
  - When DL ≤ 4, `VisualStringDistances.visual_distance(...; normalize = x -> 5 + sqrt(x))` ≤ 2.5.
- The distance check can be bypassed only with the PR label `Override AutoMerge: name similarity is okay` (L664). Registry maintainers apply that label; after that, merging is at the maintainers' discretion.

**Results.** I re-implemented these three checks verbatim with `StringDistances` 0.11 (within AutoMerge's compat `0.9, 0.10, 0.11`) and `VisualStringDistances` 0.1.1. I ran them against the 12,675 non-JLL names in General `Registry.toml` @ `0c75ff5` (2026-10-14). Script: see Appendix.

| Candidate | len | DL vs `BondGraphs` | lc-DL | visual | AutoMerge distance check |
|---|---|---|---|---|---|
| `BondGraph` | 9 | **1** | **1** | **2.27** | **FAIL** (all three checks, vs `BondGraphs`) |
| `BondGraphToolkit` | 16 | 7 | 7 | 14.15 | pass |
| `BondGraphsToolkit` | 17 | 7 | 7 | 13.92 | pass |
| `BondGraphModelingToolkit` | 24 | 15 | 15 | 29.70 | pass |
| `ModelingToolkitBondGraphs` | 25 | 15 | 15 | 43.77 | pass |
| `BondGraphDialect` | 16 | 7 | 7 | 14.05 | pass |
| `BondGraphKit` | 12 | 3 | 3 | 5.42 | pass |
| `BondGraphMTK` / `BondGraphsMTK` | 12 / 13 | 3 | 3 | 6.36 / 7.74 | pass |

What this implies:
- **`BondGraph`.** AutoMerge would block it, and a human would have to apply the override label against guideline §7, whose own example (singular vs plural of an existing name) matches this case exactly. A second, practical problem: `BondGraphs` exports a type named `BondGraph` ([src/BondGraphs.jl L46](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/BondGraphs.jl#L46)). With `using BondGraphs` and a module named `BondGraph` in the same session, the identifier is ambiguous. Our current module name is not registrable without an override, which is unlikely to be granted.
- **`BondGraphToolkit` / `BondGraphsToolkit`.** Both clear the automated checks.
  - `Toolkit` has SciML precedent for "format/domain → MTK" packages: `CellMLToolkit` ("connects CellML models to the Scientific Julia ecosystem") and `SBMLToolkit` ("…for Julia's SciML ModelingToolkit"). Their registry URLs point to SciML/CellMLToolkit.jl and SciML/SBMLToolkit.jl ([General C/CellMLToolkit/Package.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/C/CellMLToolkit/Package.toml), [S/SBMLToolkit/Package.toml](https://github.com/JuliaRegistries/General/blob/0c75ff557f058e0ce3e5734d84d350078532be3b/S/SBMLToolkit/Package.toml)). That precedent signals "MTK-based".
  - §3 slightly favors the singular `BondGraphToolkit`: we do not provide a new central type, since models are MTK systems.
  - Soft risk under §8: the Python package [BondGraphTools](https://github.com/BondGraphTools/BondGraphTools) (Apache-2.0, last push 2022-03-30, README: "BondGraphTools - A Toolkit for modelling multi-physics systems") is in the same domain. A reviewer might remark on the closeness; this is judgement, not an automated check.
  - Our own docs already call the package "BondGraphToolkit" ([docs/make.jl L139](../../../docs/make.jl): "Citing BondGraphToolkit").
- **`BondGraphModelingToolkit` / `ModelingToolkitBondGraphs`.** These pass and say "MTK-native" most explicitly, at the cost of length. `ModelingToolkit*` names in General are almost all SciML-owned: `ModelingToolkitStandardLibrary`, `…Base`, `…Tearing`, `…NeuralNets`, `…Designer`, `…Inputs`, `…Parameters`, `…Tolerances`. That prefix may read as SciML-endorsed. UNVERIFIED whether reviewers object to third-party `ModelingToolkit*` prefixes; I found no written rule.
- **`*MTK`.** These pass, and General has precedent (`ControlSystemsMTK`, `LowLevelParticleFiltersMTK`, `MTKHelpers`, `Vensim2MTK`). But `MTK` is an acronym, so it falls under §1 "avoid jargon".

**Similar registered names** (case-insensitive search of `Registry.toml`):
- Containing "bond": only `BondGraphs` and the unrelated `CrystalBonds`.
- Port-Hamiltonian neighbors: `PortHamiltonianSystems`, `PortHamiltonianBenchmarkSystems`, `PortHamiltonianModelReduction`.
- Unregistered Julia bond-graph repos on GitHub: [TRACER-LULab/SystemBondGraphs.jl](https://github.com/TRACER-LULab/SystemBondGraphs.jl) (MTK + MetaGraphsNext, no license file, last push 2024-06-28) and [mic-pan/bond_graph_model_reduction](https://github.com/mic-pan/bond_graph_model_reduction).
- Our UUID `41e83443-7082-4de0-b466-8936bcd80675` is not in General.

## Appendix: name-check script

Run with `julia --project=<env with StringDistances@0.11, VisualStringDistances@0.1> check.jl Registry.toml`:

```julia
using StringDistances, VisualStringDistances, TOML
reg = TOML.parsefile(ARGS[1])["packages"]
names = filter(!endswith("_jll"), [v["name"] for v in values(reg)])
dl(a,b) = StringDistances.DamerauLevenshtein()(a,b)
vd(a,b) = VisualStringDistances.visual_distance(a,b; normalize=x->5+sqrt(x))
for c in ["BondGraph","BondGraphToolkit","BondGraphsToolkit","BondGraphModelingToolkit"]
    hits = [n for n in names if lowercase(n) != lowercase(c) &&
            (dl(c,n) <= 2 || dl(lowercase(c),lowercase(n)) <= 1 || (dl(c,n) <= 4 && vd(c,n) <= 2.5))]
    println(c, " => ", isempty(hits) ? "PASS" : "FAIL $hits")
end
```
