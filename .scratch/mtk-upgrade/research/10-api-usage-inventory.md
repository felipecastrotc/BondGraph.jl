# Current BondGraph API and usage inventory

- Scope: the current working tree, read statically; Julia was not run.
- `src/BondGraph.jl` explicitly exports 19 names: 15 locally defined API names, three names without a local definition, and MTK's `structural_simplify` binding.
- The six maintained examples use `@named` element constructors, junctions or explicit `connect`, then `simplifysys`; generated Markdown repeats much of this code.
- Six tuple `@named` calls exist, all in two-port unit tests; examples use `coneqs` instead.
- The electrical RC example calls `Spring` for a capacitor and `Damper` for a resistor; the DC motor calls `Mass` for an inductor.
- `generate_graph` and bond expansion depend directly on MTK connection internals removed or changed in MTK 11.
- `renamevars`, `addsubsys`, and `latexify` have no local definition; the latter is called with `Symbolics.Latexify` imported by examples.

## Counting method and scope

Paths below are relative to the repo root. `E`, `T`, and `D` mean respectively `examples/*.jl`, `test/**/*.jl` (including `test/runtests.jl`), and `docs/src/**/*.md`. Each count is the number of source lines containing an active call or binding use, with its full `file:line` list; multiple calls on one line count once. In `.jl`, prose comments, test-set titles and imports are not call sites. In Markdown, fenced Julia code is a call site; a Documenter `@docs` entry is marked separately as a documentation reference. Generated docs are counted independently from their example sources, because they are separate files. The local `tmp_dev_acoustic_carol.jl` and generated `tmp_dev_*.md` are included because the requested globs include them, although `tmp_dev_*` is scratch content. This inventory covers BondGraph-owned names, including its MTK method overwrites, rather than every name re-exported wholesale by `@reexport using ModelingToolkit` (`src/BondGraph.jl:36`).

## Definitions and signatures

`src/BondGraph.jl:40-58` is the explicit export list. A single displayed signature means one local method; defaults are shown as written in source. `—` means no local definition. The `@connector` and `@oneport` entries are macros, not ordinary function methods.

| Symbol | Local signature and definition | One-line purpose |
|---|---|---|
| `Power` | `@connector function Power(; name, effort=0.0, flow=0.0, type=op)` — `src/bg.jl:27` | Make an effort/flow power port carrying bond metadata. |
| `Se` | `Se(expr; name)` — `src/bg.jl:58` | Impose a source effort. |
| `Sf` | `Sf(expr; name)` — `src/bg.jl:102` | Impose a source flow. |
| `Dq` | `Dq(; name, x=0.0)` — `src/bg.jl:128` | Make an auxiliary displacement integrator (`D(q) ~ f`); `x` is accepted but unused in the body. |
| `Junction0` | `Junction0(ps...; name="")` — `src/bg.jl:159` | Compose supplied elements into a 0-junction with directed bonds. |
| `Junction1` | `Junction1(ps...; name="")` — `src/bg.jl:207` | Compose supplied elements into a 1-junction with directed bonds. |
| `mGY` | `mGY(subsys...; name="", g=1.0, coneqs=nothing)` — `src/bg.jl:254` | Make a modulated gyrator; return `(sys, con)` when `coneqs=nothing`, otherwise append bonds and return `sys` (`src/bg.jl:284-293`). |
| `mTF` | `mTF(subsys...; name="", r=1.0, coneqs=nothing)` — `src/bg.jl:317` | Make a modulated transformer with the same two return modes (`src/bg.jl:347-356`). |
| `Mass` | `Mass(; name="", m=1.0, u=0.0)` — `src/bg.jl:379` | Make the generic inertia (I) element, with initial flow `u`. |
| `Spring` | `Spring(; name="", k=10, x=0.0)` — `src/bg.jl:406` | Make the generic compliance (C) element, with displacement state `q`. |
| `Damper` | `Damper(; name="", c=10, u=1.0)` — `src/bg.jl:434` | Make the generic resistance (R) element. |
| `isindependent` | `isindependent(var::Num)` — `src/utils/bg_utils.jl:150` | Classify a symbol as independent of time. |
| `generate_graph` | `generate_graph(mdl, var=:e; method=:stress)` — `src/utils/graph.jl:53` | Plot the bond connection graph through GraphRecipes/Plots. |
| `simplifysys` | `simplifysys(sys::ODESystem; name)` — `src/utils/system_utils.jl:16` | Expand bonds, structurally simplify, and substitute observed equations. |
| `renamevars` | —; export `src/BondGraph.jl:57`; commented old function `src/mtk.jl:232-244` | Intended variable renaming; no live local implementation. |
| `latexify` | —; export `src/BondGraph.jl:57` | Intended LaTeX rendering; examples separately import `Symbolics.Latexify`. |
| `structural_simplify` | No BondGraph method; imported from MTK at `src/utils/system_utils.jl:1`, exported at `src/BondGraph.jl:57` | MTK's structural simplification function, called as `structural_simplify(sys)`. |
| `addsubsys` | —; export `src/BondGraph.jl:57` | Intended subsystem addition; no live local implementation. |
| `@oneport` | `macro oneport(ex)` — `src/utils/macros.jl:5` | Rewrite a custom one-port function to create `power`, infer variables and compose an `ODESystem`. |

Non-exported BondGraph-owned names reached by examples or tests are below. The tag types have zero-field `struct Name end` definitions. Test imports themselves are listed later for `t`/`D`; other imports are omitted from call counts.

| Symbol | Local signature and definition | One-line purpose |
|---|---|---|
| `t` | `@variables t` — `src/bg.jl:1` | Independent time variable shared with BondGraph systems. |
| `D` | `D = Differential(t)` — `src/bg.jl:2` | Time derivative operator. |
| `expand_connections` | `expand_connections(sys::AbstractSystem, find=nothing, replace=nothing; debug=false, tol=1e-10)` — `src/mtk.jl:25-40` | MTK method overwrite that inserts bond equations during connection expansion. |
| `generate_connection_set` | `generate_connection_set(sys::AbstractSystem, find=nothing, replace=nothing)` — `src/mtk.jl:42-54` | MTK method overwrite that extracts bond connection sets before merging. |
| `substitute_dict` | `substitute_dict(O, expr, sub)` — `src/utils/system_utils.jl:61` | Recursively substitute keys found in a symbolic expression. |
| `set_bg_metadata` | `set_bg_metadata(s, type)` — `src/types.jl:18` | Attach bond type and role metadata to a symbolic variable. |
| `get_bg_junction` | `get_bg_junction(s)` — `src/types.jl:19` | Read that bond metadata. |
| `update_mtk_con` | `update_mtk_con(s, type)` — `src/types.jl:20` | Set MTK `VariableConnectType` metadata on a symbol. |
| `flatinput` | `flatinput(ps::AbstractArray)` — `src/utils/bg_utils.jl:96` | Flatten junction arguments into systems and sign values. |
| `equalityeqs` | `equalityeqs(con::AbstractArray)` — `src/utils/bg_utils.jl:26` | Generate equalities among port variables. |
| `sumvar` | `sumvar(con::AbstractArray)` — `src/utils/bg_utils.jl:67` | Generate a zero-sum junction equation. |
| `check_bg_con` | `check_bg_con(connectionset::ConnectionSet)` — `src/utils/connection_utils.jl:18` | Check whether an MTK connection set is a bond. |
| `get_bg_connection_set!` | `get_bg_connection_set!(connectionsets::Array{ConnectionSet})` — `src/utils/connection_utils.jl:35` | Remove and return bond sets from MTK's connection sets. |
| `add_idx` | `add_idx(var, idx)` — `src/utils/connection_utils.jl:61` | Rename a symbolic variable with a numeric suffix. |
| `gen_tp_con!` | `gen_tp_con!(eqs::AbstractArray, sys::ODESystem, subsys)` — `src/utils/connection_utils.jl:89` | Append bonds for a two-port's input/output systems. |
| `csets2dict` | `csets2dict(csets)` — `src/utils/graph.jl:109` | Index connection elements by namespaced variable name. |
| `csets2adjmtx` | `csets2adjmtx(csets, str2con; filterstr="f(t)", filterflow=false)` — `src/utils/graph.jl:144` | Turn directed connection sets into an adjacency matrix. |
| `bg`, `bgeffort`, `bgflow` | `struct … end` — `src/types.jl:7-9` | Metadata tags for bond connection and effort/flow roles. |
| `op`, `tpgy`, `j0`, `j1` | `struct … end` — `src/types.jl:10-11,13-14` | Metadata tags for one-port, gyrator, and junction kinds. |

## Call sites

The following ledger gives **all** active call sites for the explicit exports and for non-exported functions used by examples/tests. `D` includes fenced code and `@docs` entries; the latter are marked `[ref]`. A `0` means no call in the requested globs.

### `Power`

- E (0): —
- T (29): `test/unit-tests/test_bg_utils.jl:90`; `test/unit-tests/test_connection_utils.jl:20,68,69`; `test/unit-tests/test_graph.jl:19,20`; `test/unit-tests/test_macros.jl:33`; `test/unit-tests/test_ports.jl:15,20,24,30,36,59,81,110,132,161,189,218,255,312,377,378,394,395,431,432,448,449`
- D (1): `docs/src/lib/public.md:39[ref]`

### `Se`

- E (6): `examples/electrical_dc_motor.jl:66,71`; `examples/mechanical_duffing.jl:94`; `examples/mechanical_msd_1dof.jl:56`; `examples/mechanical_msd_2dof.jl:60`; `examples/tmp_dev_acoustic_carol.jl:34`
- T (9): `test/regression/test_electrical_dc_motor.jl:6,9`; `test/regression/test_mechanical_duffing.jl:18`; `test/regression/test_mechanical_msd_1dof.jl:7`; `test/regression/test_mechanical_msd_2dof.jl:7`; `test/unit-tests/test_ports.jl:53,56,75,78`
- D (8): `docs/src/examples/electrical_dc_motor.md:76,84`; `docs/src/examples/mechanical_duffing.md:112`; `docs/src/examples/mechanical_msd_1dof.md:70`; `docs/src/examples/mechanical_msd_2dof.md:74`; `docs/src/examples/tmp_dev_acoustic_carol.md:39`; `docs/src/examples/tmp_dev_chat.md:59`; `docs/src/lib/public.md:30[ref]`

### `Sf`

- E (0): —
- T (4): `test/unit-tests/test_ports.jl:104,107,126,129`
- D (1): `docs/src/lib/public.md:31[ref]`

### `Dq`

- E (0): —
- T (0): —
- D (0): —

### `Junction0`

- E (6): `examples/electrical_rc.jl:36`; `examples/mechanical_msd_2dof.jl:71`; `examples/mechanical_msd_ndof.jl:58`; `examples/tmp_dev_acoustic_carol.jl:50,56,62`
- T (5): `test/regression/test_electrical_rc.jl:5`; `test/regression/test_mechanical_msd_2dof.jl:10`; `test/regression/test_mechanical_msd_ndof.jl:10`; `test/unit-tests/test_ports.jl:246,249`
- D (8): `docs/src/examples/mechanical_msd_2dof.md:91`; `docs/src/examples/mechanical_msd_ndof.md:81`; `docs/src/examples/tmp_dev_acoustic_carol.md:55,61,67`; `docs/src/getting_started_auto.md:61,132`; `docs/src/lib/public.md:33[ref]`

### `Junction1`

- E (13): `examples/electrical_dc_motor.jl:77,78`; `examples/mechanical_duffing.jl:99`; `examples/mechanical_msd_1dof.jl:61`; `examples/mechanical_msd_2dof.jl:67,70,72`; `examples/mechanical_msd_ndof.jl:47,50,59`; `examples/tmp_dev_acoustic_carol.jl:49,53,59`
- T (13): `test/regression/test_electrical_dc_motor.jl:10,11`; `test/regression/test_mechanical_duffing.jl:19`; `test/regression/test_mechanical_msd_1dof.jl:8`; `test/regression/test_mechanical_msd_2dof.jl:8,9,11`; `test/regression/test_mechanical_msd_ndof.jl:6,7,11`; `test/unit-tests/test_ports.jl:303,306`; `test/unit-tests/test_system_utils.jl:15`
- D (15): `docs/src/examples/electrical_dc_motor.md:96,97`; `docs/src/examples/mechanical_duffing.md:120`; `docs/src/examples/mechanical_msd_1dof.md:78`; `docs/src/examples/mechanical_msd_2dof.md:84,90,92`; `docs/src/examples/mechanical_msd_ndof.md:60,66,82`; `docs/src/examples/tmp_dev_acoustic_carol.md:54,58,64`; `docs/src/examples/tmp_dev_chat.md:66`; `docs/src/lib/public.md:32[ref]`

### `mGY`

- E (1): `examples/electrical_dc_motor.jl:82`
- T (5): `test/regression/test_electrical_dc_motor.jl:13`; `test/unit-tests/test_ports.jl:361,367,372,389`
- D (2): `docs/src/examples/electrical_dc_motor.md:104`; `docs/src/lib/public.md:34[ref]`

### `mTF`

- E (1): `examples/tmp_dev_acoustic_carol.jl:63`
- T (4): `test/unit-tests/test_ports.jl:415,421,426,443`
- D (2): `docs/src/examples/tmp_dev_acoustic_carol.md:68`; `docs/src/lib/public.md:35[ref]`

### `Mass`

- E (9): `examples/electrical_dc_motor.jl:64,69`; `examples/mechanical_duffing.jl:87`; `examples/mechanical_msd_1dof.jl:50`; `examples/mechanical_msd_2dof.jl:54`; `examples/mechanical_msd_ndof.jl:38`; `examples/tmp_dev_acoustic_carol.jl:38,43,46`
- T (16): `test/regression/test_electrical_dc_motor.jl:4,7`; `test/regression/test_mechanical_duffing.jl:13`; `test/regression/test_mechanical_msd_1dof.jl:3`; `test/regression/test_mechanical_msd_2dof.jl:3`; `test/regression/test_mechanical_msd_ndof.jl:3`; `test/unit-tests/test_bg_utils.jl:63`; `test/unit-tests/test_connection_utils.jl:16,87`; `test/unit-tests/test_ports.jl:156,158,243,300,357,411`; `test/unit-tests/test_system_utils.jl:13`
- D (11): `docs/src/examples/electrical_dc_motor.md:74,82`; `docs/src/examples/mechanical_duffing.md:102`; `docs/src/examples/mechanical_msd_1dof.md:60`; `docs/src/examples/mechanical_msd_2dof.md:64`; `docs/src/examples/mechanical_msd_ndof.md:48`; `docs/src/examples/tmp_dev_acoustic_carol.md:43,48,51`; `docs/src/examples/tmp_dev_chat.md:50`; `docs/src/lib/public.md:36[ref]`

### `Spring`

- E (7): `examples/electrical_rc.jl:32`; `examples/mechanical_duffing.jl:88`; `examples/mechanical_msd_1dof.jl:51`; `examples/mechanical_msd_2dof.jl:55`; `examples/mechanical_msd_ndof.jl:39`; `examples/tmp_dev_acoustic_carol.jl:40,41`
- T (7): `test/regression/test_electrical_rc.jl:3`; `test/regression/test_mechanical_duffing.jl:14`; `test/regression/test_mechanical_msd_1dof.jl:4`; `test/regression/test_mechanical_msd_2dof.jl:4`; `test/regression/test_mechanical_msd_ndof.jl:4`; `test/unit-tests/test_ports.jl:184,186`
- D (10): `docs/src/examples/mechanical_duffing.md:103`; `docs/src/examples/mechanical_msd_1dof.md:61`; `docs/src/examples/mechanical_msd_2dof.md:65`; `docs/src/examples/mechanical_msd_ndof.md:49`; `docs/src/examples/tmp_dev_acoustic_carol.md:45,46`; `docs/src/examples/tmp_dev_chat.md:51`; `docs/src/getting_started_auto.md:57,123`; `docs/src/lib/public.md:37[ref]`

### `Damper`

- E (9): `examples/electrical_dc_motor.jl:65,70`; `examples/electrical_rc.jl:33`; `examples/mechanical_duffing.jl:90`; `examples/mechanical_msd_1dof.jl:52`; `examples/mechanical_msd_2dof.jl:56`; `examples/mechanical_msd_ndof.jl:40`; `examples/tmp_dev_acoustic_carol.jl:35,36`
- T (16): `test/regression/test_electrical_dc_motor.jl:5,8`; `test/regression/test_electrical_rc.jl:4`; `test/regression/test_mechanical_duffing.jl:16`; `test/regression/test_mechanical_msd_1dof.jl:5`; `test/regression/test_mechanical_msd_2dof.jl:5`; `test/regression/test_mechanical_msd_ndof.jl:5`; `test/unit-tests/test_bg_utils.jl:64`; `test/unit-tests/test_connection_utils.jl:17`; `test/unit-tests/test_ports.jl:213,215,244,301,358,412`; `test/unit-tests/test_system_utils.jl:14`
- D (12): `docs/src/examples/electrical_dc_motor.md:75,83`; `docs/src/examples/mechanical_duffing.md:105`; `docs/src/examples/mechanical_msd_1dof.md:62`; `docs/src/examples/mechanical_msd_2dof.md:66`; `docs/src/examples/mechanical_msd_ndof.md:50`; `docs/src/examples/tmp_dev_acoustic_carol.md:40,41`; `docs/src/examples/tmp_dev_chat.md:52`; `docs/src/getting_started_auto.md:58,124`; `docs/src/lib/public.md:38[ref]`

### `isindependent`

- E (0): —
- T (4): `test/unit-tests/test_bg_utils.jl:92,93,95,96`
- D (1): `docs/src/lib/public.md:47[ref]`

### `generate_graph`

- E (7): `examples/electrical_dc_motor.jl:95`; `examples/electrical_rc.jl:54`; `examples/mechanical_duffing.jl:104`; `examples/mechanical_msd_1dof.jl:65`; `examples/mechanical_msd_2dof.jl:95`; `examples/mechanical_msd_ndof.jl:85`; `examples/tmp_dev_acoustic_carol.jl:70`
- T (0): —
- D (9): `docs/src/examples/electrical_dc_motor.md:126`; `docs/src/examples/mechanical_duffing.md:128`; `docs/src/examples/mechanical_msd_1dof.md:86`; `docs/src/examples/mechanical_msd_2dof.md:128`; `docs/src/examples/mechanical_msd_ndof.md:128`; `docs/src/examples/tmp_dev_acoustic_carol.md:75`; `docs/src/examples/tmp_dev_chat.md:73`; `docs/src/getting_started_auto.md:81`; `docs/src/lib/public.md:45[ref]`

### `simplifysys`

- E (7): `examples/electrical_dc_motor.jl:101`; `examples/electrical_rc.jl:38`; `examples/mechanical_duffing.jl:111`; `examples/mechanical_msd_1dof.jl:71`; `examples/mechanical_msd_2dof.jl:101`; `examples/mechanical_msd_ndof.jl:91`; `examples/tmp_dev_acoustic_carol.jl:74`
- T (7): `test/regression/test_electrical_dc_motor.jl:16`; `test/regression/test_electrical_rc.jl:6`; `test/regression/test_mechanical_duffing.jl:20`; `test/regression/test_mechanical_msd_1dof.jl:9`; `test/regression/test_mechanical_msd_2dof.jl:19`; `test/regression/test_mechanical_msd_ndof.jl:19`; `test/unit-tests/test_system_utils.jl:18`
- D (10): `docs/src/examples/electrical_dc_motor.md:138`; `docs/src/examples/mechanical_duffing.md:141`; `docs/src/examples/mechanical_msd_1dof.md:98`; `docs/src/examples/mechanical_msd_2dof.md:140`; `docs/src/examples/mechanical_msd_ndof.md:140`; `docs/src/examples/tmp_dev_acoustic_carol.md:79`; `docs/src/examples/tmp_dev_chat.md:85`; `docs/src/getting_started_auto.md:63,133`; `docs/src/lib/public.md:46[ref]`

### `renamevars`

- E (0): —
- T (0): —
- D (0): —

### `latexify`

- E (6): `examples/electrical_dc_motor.jl:115`; `examples/mechanical_duffing.jl:108,125`; `examples/mechanical_msd_1dof.jl:85`; `examples/mechanical_msd_2dof.jl:115`; `examples/tmp_dev_acoustic_carol.jl:82`
- T (0): —
- D (7): `docs/src/examples/electrical_dc_motor.md:164`; `docs/src/examples/mechanical_duffing.md:135,168`; `docs/src/examples/mechanical_msd_1dof.md:124`; `docs/src/examples/mechanical_msd_2dof.md:166`; `docs/src/examples/tmp_dev_acoustic_carol.md:87`; `docs/src/examples/tmp_dev_chat.md:99`

### `structural_simplify`

- E (1): `examples/electrical_dc_motor.jl:159`
- T (0): —
- D (1): `docs/src/examples/electrical_dc_motor.md:230`

### `addsubsys`

- E (0): —
- T (0): —
- D (0): —

### `@oneport`

- E (1): `examples/mechanical_duffing.jl:55`
- T (2): `test/regression/test_mechanical_duffing.jl:2`; `test/unit-tests/test_macros.jl:12`
- D (1): `docs/src/examples/mechanical_duffing.md:64`

### `expand_connections`

- E (7): `examples/electrical_dc_motor.jl:98`; `examples/mechanical_duffing.jl:107,108`; `examples/mechanical_msd_1dof.jl:68`; `examples/mechanical_msd_2dof.jl:98`; `examples/mechanical_msd_ndof.jl:88`; `examples/tmp_dev_acoustic_carol.jl:72`
- T (2): `test/unit-tests/test_ports.jl:285,342`
- D (8): `docs/src/examples/electrical_dc_motor.md:132`; `docs/src/examples/mechanical_duffing.md:134,135`; `docs/src/examples/mechanical_msd_1dof.md:92`; `docs/src/examples/mechanical_msd_2dof.md:134`; `docs/src/examples/mechanical_msd_ndof.md:134`; `docs/src/examples/tmp_dev_acoustic_carol.md:77`; `docs/src/examples/tmp_dev_chat.md:79`

### `generate_connection_set`

- E (0): —
- T (1): `test/unit-tests/test_graph.jl:31`
- D (0): —

### `substitute_dict`

- E (0): —
- T (2): `test/unit-tests/test_system_utils.jl:38,42`
- D (0): —

### `set_bg_metadata`

- E (0): —
- T (2): `test/unit-tests/test_types.jl:12,18`
- D (0): —

### `get_bg_junction`

- E (0): —
- T (7): `test/unit-tests/test_ports.jl:39,43,252,253,309,310`; `test/unit-tests/test_types.jl:19`
- D (0): —

### `update_mtk_con`

- E (0): —
- T (1): `test/unit-tests/test_types.jl:25`
- D (0): —

### `flatinput`

- E (0): —
- T (3): `test/unit-tests/test_bg_utils.jl:51,67,82`
- D (0): —

### `equalityeqs`

- E (0): —
- T (3): `test/unit-tests/test_bg_utils.jl:16,30,45`
- D (0): —

### `sumvar`

- E (0): —
- T (1): `test/unit-tests/test_bg_utils.jl:40`
- D (0): —

### `check_bg_con`

- E (0): —
- T (1): `test/unit-tests/test_connection_utils.jl:37`
- D (0): —

### `get_bg_connection_set!`

- E (0): —
- T (1): `test/unit-tests/test_connection_utils.jl:41`
- D (0): —

### `add_idx`

- E (0): —
- T (2): `test/unit-tests/test_connection_utils.jl:52,53`
- D (0): —

### `gen_tp_con!`

- E (0): —
- T (6): `test/unit-tests/test_connection_utils.jl:76,77,82,92,98,104`
- D (0): —

### `csets2dict`

- E (0): —
- T (6): `test/unit-tests/test_graph.jl:36,41,54,61,69,81`
- D (0): —

### `csets2adjmtx`

- E (0): —
- T (4): `test/unit-tests/test_graph.jl:64,71,83,93`
- D (0): —

## Non-exported bindings and metadata tags: use sites

`t` and `D` are not in the explicit export list, but are defined in `src/bg.jl:1-2`; callers import them explicitly. This ledger omits imports, comments, string literals and property access such as `sol.t`. It includes bindings inside symbolic declarations (`x(t)`) and `D(x)` calls.

| Binding | E | T | D |
|---|---|---|---|
| `t` | 17: `examples/electrical_dc_motor.jl:87,134,135,136,156`; `examples/mechanical_duffing.jl:58,93,146`; `examples/mechanical_msd_1dof.jl:55,105`; `examples/mechanical_msd_2dof.jl:59,88,135`; `examples/mechanical_msd_ndof.jl:79`; `examples/tmp_dev_acoustic_carol.jl:9,66,79` | 31: `test/regression/test_electrical_dc_motor.jl:14`; `test/regression/test_mechanical_duffing.jl:3,17,22`; `test/regression/test_mechanical_msd_1dof.jl:6,11`; `test/regression/test_mechanical_msd_2dof.jl:6,17,21`; `test/regression/test_mechanical_msd_ndof.jl:17`; `test/unit-tests/test_bg_utils.jl:12,36,87,88`; `test/unit-tests/test_connection_utils.jl:25,48,71`; `test/unit-tests/test_graph.jl:27`; `test/unit-tests/test_macros.jl:15,31`; `test/unit-tests/test_ports.jl:73,124,190,273,274,330,331`; `test/unit-tests/test_system_utils.jl:33`; `test/unit-tests/test_types.jl:11,17,24` | 19: `docs/src/examples/electrical_dc_motor.md:112,194,195,196,224`; `docs/src/examples/mechanical_duffing.md:67,111,204`; `docs/src/examples/mechanical_msd_1dof.md:69,159`; `docs/src/examples/mechanical_msd_2dof.md:73,114,201`; `docs/src/examples/mechanical_msd_ndof.md:115`; `docs/src/examples/tmp_dev_acoustic_carol.md:14,71,84`; `docs/src/examples/tmp_dev_chat.md:58,126` |
| `D` | 4: `examples/electrical_dc_motor.jl:150,151,152`; `examples/mechanical_duffing.jl:65` | 6: `test/regression/test_mechanical_duffing.jl:7`; `test/unit-tests/test_macros.jl:22,38`; `test/unit-tests/test_ports.jl:168,197`; `test/unit-tests/test_system_utils.jl:21` | 4: `docs/src/examples/electrical_dc_motor.md:215,216,217`; `docs/src/examples/mechanical_duffing.md:74` |

Explicit `using`/`import BondGraph: t, D` sites: `examples/electrical_dc_motor.jl:55`, `examples/mechanical_duffing.jl:43`, `examples/mechanical_msd_1dof.jl:41`, `examples/mechanical_msd_2dof.jl:45`, `examples/mechanical_msd_ndof.jl:27`, `examples/tmp_dev_acoustic_carol.jl:2`, `test/regression/helpers.jl:3`, `test/unit-tests/test_graph.jl:2`, `test/unit-tests/test_macros.jl:2`, `test/unit-tests/test_ports.jl:2`, `test/unit-tests/test_system_utils.jl:6`, `test/unit-tests/test_types.jl:2`. `test/unit-tests/test_bg_utils.jl:6` and `test/unit-tests/test_connection_utils.jl:6` import `t` alone. Generated Markdown duplicates the six example imports at `docs/src/examples/electrical_dc_motor.md:62`, `mechanical_duffing.md:50`, `mechanical_msd_1dof.md:48`, `mechanical_msd_2dof.md:52`, `mechanical_msd_ndof.md:34`, `tmp_dev_acoustic_carol.md:7` (all under `docs/src/examples/`).

Tag-type uses in tests (imports excluded): `bg` — `test/unit-tests/test_types.jl:13`; `bgeffort` — `test/unit-tests/test_ports.jl:41`; `bgflow` — `test/unit-tests/test_ports.jl:45`; `op` — `test/unit-tests/test_types.jl:12,13,25,26`; `tpgy` — `test/unit-tests/test_connection_utils.jl:68,69`; `j0` — `test/unit-tests/test_ports.jl:36,40,44,252,253,255`, `test/unit-tests/test_types.jl:18,20`; `j1` — `test/unit-tests/test_connection_utils.jl:20`, `test/unit-tests/test_ports.jl:309,310,312`. The `j0`/`j1` names in N-DoF examples and regression tests are local junction variables, not the metadata tag types. `tptf` (`src/types.jl:12`) is used internally by `mTF` but not directly by examples/tests.

## How models are built today

These three snippets are verbatim from the example source files. They cover a single junction, explicit junction-to-junction bonds, and a custom nonlinear compliance.

**RC circuit — `examples/electrical_rc.jl:31-46`:**

```julia
## Define the elements of the RC circuit
@named C = Spring(k=1 / 1.0)   # Electric capacitance (F)
@named R = Damper(c=1.0)       # Electric resistance (Ohm)

## Create the RC circuit model by connecting the elements
@named model = Junction0(R, C)
## Simplify the system equations
@named sys = simplifysys(model)

## Set the time span for simulation
tspan = (0, 10)
## Define the initial conditions
initial_conditions = [R.power.f => 0.0, C.q => 1.0]

## Define the ODE problem
prob = ODEProblem(sys, initial_conditions, tspan)
```

**Two degrees of freedom — `examples/mechanical_msd_2dof.jl:67-90`:**

```julia
@named m1_j1 = Junction1([-1, m], [-1, d], [-1, s])

# Then, we define the $2^{nd}$ DoF mass (($m_2$)) with two 1-junctions and one 0-junction.
@named m2_j1 = Junction1([-1, m], f)
@named m2_j0 = Junction0()
@named m2_sd = Junction1([-1, s], [-1, d])

# ## Model Construction
# 
# ### Connect the subsystems

# Define the connections between the junctions.
cons = [
    connect(m2_j0.power, m2_j1.power),
    connect(m2_j0.power, m2_sd.power),
    connect(m1_j1.power, m2_j0.power),
]

# ### Build the system 
# 
# We build the system by first creating a `ODESystem` with junctions connections and setting the independent variable as t.
@named mdl = ODESystem(cons, t)
# Then, we need to add the junctions subsystems to the model
mdl = compose(mdl, m1_j1, m2_j0, m2_j1, m2_sd)
```

**Custom cubic compliance — `examples/mechanical_duffing.jl:55-67`:**

```julia
@oneport function Spring3(; name, k = 1.0, x = 0.0)

    ## Initialize the displacement state with its initial value `x`
    @variables q(t) = x
    ## Set the equation parameters
    @parameters C = 1 / k

    ## Define Spring3 element equations with the non-linearity of the cubic spring
    [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]
end
```

Other concrete patterns:

- **Two-port return modes:** `examples/electrical_dc_motor.jl:81-90` builds `eqs = []`, then `@named gy = mGY(je, jm, g = g, coneqs = eqs)`, then `ODESystem(eqs, t)` and `compose(mdl, gy, je, jm)`. The regression copy is `test/regression/test_electrical_dc_motor.jl:13-16`; generated docs are `docs/src/examples/electrical_dc_motor.md:103-117`. The scratch acoustic example uses `@named tf = mTF(j01, j03, r=2, coneqs=cons)` at `examples/tmp_dev_acoustic_carol.jl:63` (generated Markdown `docs/src/examples/tmp_dev_acoustic_carol.md:68`). In contrast, the default `coneqs=nothing` returns `(sys, con)`.
- **Every tuple `@named` use:** `test/unit-tests/test_ports.jl:361,372,389` call `@named sys, con/eqs = mGY(...)`; `test/unit-tests/test_ports.jl:415,426,443` do the same for `mTF`. The same file also uses ordinary tuple assignment with explicit `name` at `:367` and `:421`. No tuple `@named` call occurs in `examples/*.jl` or `docs/src/**/*.md`. The tuple form relies on the `split_assign` and `single_named_expr` MTK method overwrites in `src/mtk.jl:59-92`, as the research note explains.
- **Custom element:** `@oneport` appears at `examples/mechanical_duffing.jl:55`, `test/regression/test_mechanical_duffing.jl:2`, `test/unit-tests/test_macros.jl:12`, and `docs/src/examples/mechanical_duffing.md:64`; the macro supplies a `power` port and infers symbolic states/parameters (`src/utils/macros.jl:5-42`). The same `Spring3` is then constructed and passed to a junction at `examples/mechanical_duffing.jl:89,99` and the regression test's `:15,19`. `Power()` itself is used directly in tests, never called directly in an example.
- **Signs:** there is no `-x` unary-system argument. Direction is requested with two-element vectors such as `[-1, m]`, parsed by `flatinput` (`src/utils/bg_utils.jl:96-135`); an unsigned system means `+1`. Every signed junction construction in examples is at `examples/electrical_dc_motor.jl:77,78`, `mechanical_duffing.jl:99`, `mechanical_msd_1dof.jl:61`, `mechanical_msd_2dof.jl:67,70,72`, `mechanical_msd_ndof.jl:47,50,59`, `tmp_dev_acoustic_carol.jl:49,53,56,59,62` (all under `examples/`). Regression copies are `test/regression/test_electrical_dc_motor.jl:10,11`, `test_mechanical_duffing.jl:19`, `test_mechanical_msd_1dof.jl:8`, `test_mechanical_msd_2dof.jl:8,9,11`, `test_mechanical_msd_ndof.jl:6,7,11` (all under `test/regression/`). Markdown copies are `docs/src/examples/electrical_dc_motor.md:96,97`, `mechanical_duffing.md:120`, `mechanical_msd_1dof.md:78`, `mechanical_msd_2dof.md:84,90,92`, `mechanical_msd_ndof.md:60,66,82`, `tmp_dev_acoustic_carol.md:54,58,61,64,67`, `tmp_dev_chat.md:66` (all under `docs/src/examples/`). `test/unit-tests/test_bg_utils.jl:76-82` also tests arbitrary signed numeric prefixes including `[-10, mass]` and `[-1.0, damper]` through `flatinput`.
- **`couple`:** zero active keyword uses in the requested files. The signatures at `src/bg.jl:159,207` accept only `name`; only their docstring signatures at `src/bg.jl:138,186` claim `couple=true`. The two-port docstrings also claim `coneqs=[]` (`src/bg.jl:234,297`), whereas the implementation defaults to `nothing` (`src/bg.jl:254,317`).
- **Compilation, visualization and observation:** all six maintained examples call `simplifysys`; five also show `expand_connections` and all six call `generate_graph` (the RC source comments on the graph at `examples/electrical_rc.jl:52` and calls it at `:54`). `simplifysys` calls `reducedobs` internally (`src/utils/system_utils.jl:16-20`); no example or active test calls `reducedobs` directly. `substitute_dict` has two active unit-test calls (`test/unit-tests/test_system_utils.jl:38,42`) and is called by `substitute_observed` (`src/utils/system_utils.jl:105,110`), which `reducedobs` calls. The observed-substitution test is only a TODO comment at `test/unit-tests/test_system_utils.jl:45-47`.
- **Time and derivatives:** `using BondGraph: t, D` is repeated in the five non-RC maintained examples (`examples/electrical_dc_motor.jl:55`, `mechanical_duffing.jl:43`, `mechanical_msd_1dof.jl:41`, `mechanical_msd_2dof.jl:45`, `mechanical_msd_ndof.jl:27`); the acoustic scratch example does so at `:2`. Tests import these bindings too, as listed above. Examples use `t` in `@parameters F(t)`, `@variables x(t)` and `ODESystem(..., t)`, and use `D(q)` in custom equations. The `@parameters F(t)` sites are `examples/mechanical_duffing.jl:93`, `mechanical_msd_1dof.jl:55`, `mechanical_msd_2dof.jl:59`, `tmp_dev_acoustic_carol.jl:9`; the research note says this syntax is an error in MTK 11.

## Mechanical names for generic I/C/R elements

The constructor bodies already store the domain-neutral parameters `I`, `C`, and `R` (`src/bg.jl:381,410,437`). Their public names are mechanical: `Mass` = inertia, `Spring` = compliance, and `Damper` = resistance. The complete constructor-site ledger above lists every code occurrence across the requested globs, including mechanical uses. The cross-domain uses that expose the naming mismatch are:

| Domain and physical role | Source uses | Test / generated-doc copies |
|---|---|---|
| Electrical capacitor as `Spring` | `examples/electrical_rc.jl:32` (`# Electric capacitance (F)`) | `test/regression/test_electrical_rc.jl:3`; `docs/src/getting_started_auto.md:57,123` |
| Electrical resistor as `Damper` | `examples/electrical_rc.jl:33`; `examples/electrical_dc_motor.jl:65` (armature resistance) | `test/regression/test_electrical_rc.jl:4`; `test/regression/test_electrical_dc_motor.jl:5`; `docs/src/getting_started_auto.md:58,124`; `docs/src/examples/electrical_dc_motor.md:75` |
| Electrical inductor as `Mass` | `examples/electrical_dc_motor.jl:64` (`# (H) Electric inductance`) | `test/regression/test_electrical_dc_motor.jl:4`; `docs/src/examples/electrical_dc_motor.md:74` |
| Mechanical rotor inertia / friction in the same electrical example | `examples/electrical_dc_motor.jl:69-70` | `test/regression/test_electrical_dc_motor.jl:7-8`; `docs/src/examples/electrical_dc_motor.md:82-83` |
| Acoustic impedance/inertia/compliance in scratch example | `examples/tmp_dev_acoustic_carol.jl:35-36` (`Damper`), `:38,43,46` (`Mass`), `:40-41` (`Spring`) | `docs/src/examples/tmp_dev_acoustic_carol.md:40-41,43,45-46,48,51` |

`docs/src/examples/tmp_dev_chat.md:50-52` repeats the mechanical 1-DoF constructors in a scratch draft. The public docs list only the mechanical names at `docs/src/lib/public.md:36-38`; the domain-neutral `I`/`C`/`R` names from `GLOSSARY.md` are not constructors in this module.

## Dead, dangling, and test-only surface

- `renamevars` and `addsubsys` are exported at `src/BondGraph.jl:57`, with no live local definition and zero requested call sites. `renamevars` appears only as commented old code at `src/mtk.jl:232-244`. `latexify` is likewise exported without a local definition; six example calls and seven generated-doc calls are recorded above. Those examples import `Symbolics.Latexify` separately (`examples/electrical_dc_motor.jl:59`, `mechanical_duffing.jl:47`, `mechanical_msd_1dof.jl:45`, `mechanical_msd_2dof.jl:49`, `mechanical_msd_ndof.jl:31`, `tmp_dev_acoustic_carol.jl:6`). This is a static binding inventory, not a runtime claim that `BondGraph.latexify` resolves.
- `structural_simplify` is a different case: it is imported from MTK, not locally defined (`src/utils/system_utils.jl:1`), and is called directly at `examples/electrical_dc_motor.jl:159` and indirectly inside `simplifysys` (`src/utils/system_utils.jl:18`). The research note says MTK 11 retains it as a deprecated alias of `mtkcompile`.
- `Dq` is exported but has no example, test or local-doc call. `GenericDamper` (`src/bg.jl:460`) is non-exported and has zero example/test calls, despite its Documenter entry at `docs/src/lib/public.md:48`; its own `src/bg.jl:461` TODO questions whether it works. `It` (`src/bg.jl:3`) is a non-exported, unused `Integral` binding. `setinitialval` (`src/utils/macros.jl:46`) has no active call in `src/`, examples or tests. `update_mtk_con` and `add_idx` have test calls but no production call in included source. `tptf` is internal to `mTF`; no caller constructs the tag directly.
- The graph helpers `adjmtx2eqs`, `retrieve_node_details`, and `gen_con_name` are used in `src/utils/graph.jl` (`:32,207,238,241` as applicable), but their names occur only in TODO comments at `test/unit-tests/test_graph.jl:100-102`, not active tests. The `ConnectionSet` and `ConnectionElement` imports used by that file (`:6-7`) and `src/utils/graph.jl:4-5` are removed in MTK 11 per the research note.

### TODO markers in `src/`

| Location | Marker |
|---|---|
| `src/bg.jl:44,88` | Source-element docstrings say `TODO: Improve` / `TODO: Improvess` (docstring text, included for completeness). |
| `src/bg.jl:255-256,318-319` | Two-port compose ordering and constitutive-sign doubts. |
| `src/bg.jl:461` | Check whether `GenericDamper` works. |
| `src/mtk.jl:95` | Check legacy commented functions. |
| `src/utils/bg_utils.jl:27,68,97` | Restrict the accepted variable/array type. |
| `src/utils/graph.jl:75` | Replace index iteration based on `length`/`size`. |
| `src/utils/macros.jl:4,44` | Add initial-state support; question `setinitialval`. |
| `src/mbg.jl:3` | Multibond experiment is old and not working; this file is not included by `src/BondGraph.jl:69`. |

## Export decision matrix

“Documented” means a direct `@docs` entry in `docs/src/lib/public.md`, not a passing example or an MTK manual entry. “MTK internal affected” is based on `.scratch/mtk-upgrade/research-mtk-internals.md`: **direct** means the symbol's body calls or imports a changed/removed MTK internal; **via bonds** means it constructs a `Power`/bond model whose current expansion depends on the overwritten connection methods. A public MTK deprecation is distinguished from an internal dependency.

| Symbol | Used by examples? | Used by tests? | Documented? | MTK internal affected in 11? |
|---|---|---|---|---|
| `Power` | No direct call | Yes | Yes `:39` | Yes, via custom `connect=bg` connection semantics and `VariableConnectType`; MTK connection representation changed. |
| `Se` | Yes | Yes | Yes `:30` | Yes, via bonds; also calls internal `isvariable`. |
| `Sf` | No | Yes | Yes `:31` | Yes, via bonds; also calls internal `isvariable`. |
| `Dq` | No | No | No | No changed connection internal; its `ODESystem` spelling is a deprecated public alias. |
| `Junction0` | Yes | Yes | Yes `:33` | Yes, via bond connection expansion. |
| `Junction1` | Yes | Yes | Yes `:32` | Yes, via bond connection expansion. |
| `mGY` | Yes | Yes | Yes `:34` | Yes, via bonds; tuple `@named` depends on MTK method overwrites. |
| `mTF` | Yes, scratch example | Yes | Yes `:35` | Yes, via bonds; tuple `@named` depends on MTK method overwrites. |
| `Mass` | Yes | Yes | Yes `:36` | Yes, via `Power`/bonds. |
| `Spring` | Yes | Yes | Yes `:37` | Yes, via `Power`/bonds. |
| `Damper` | Yes | Yes | Yes `:38` | Yes, via `Power`/bonds. |
| `isindependent` | No | Yes | Yes `:47` | Yes: internal `isvariable`, plus `istree` imported from MTK although removed there. |
| `generate_graph` | Yes | No active test | Yes `:45` | Yes, directly: `ConnectionSet`, `ConnectionElement`, `generate_connection_set!`, and `namespaced_var`. |
| `simplifysys` | Yes | Yes | Yes `:46` | Yes, directly: overwritten `expand_connections` and old `reducedobs` path. |
| `renamevars` | No | No | No | N/A: no live definition. |
| `latexify` | Yes, via separate import | No | No | N/A: no local definition; Symbolics 7 makes Latexify a weak dependency. |
| `structural_simplify` | Yes | No active call | No | No BondGraph MTK internal; MTK's public function is a deprecated alias of `mtkcompile`. |
| `addsubsys` | No | No | No | N/A: no live definition. |
| `@oneport` | Yes | Yes | No | Yes, via generated `Power`/bonds and `isindependent`. |

The research note further identifies the direct compatibility breaks behind those rows: `ConnectionSet`, `ConnectionElement` and `namespaced_var` removed in MTK 10.4; `expand_connections`/`generate_connection_set` signatures changed; Julia ≥1.10 rejects the MTK method overwrites during precompilation; `ODESystem` and `structural_simplify` remain deprecated aliases; and `reducedobs` constructs an `ODESystem(eqs; name)` form that the note says no longer works. This is an inventory for the API decision ticket, not a migration verdict.

STATUS: DONE
