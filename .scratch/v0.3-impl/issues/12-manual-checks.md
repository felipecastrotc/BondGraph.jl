# 12: Manual checks before the release

**What to build:** Nothing to code. The user checks these by hand, then decides whether anything changes.

**Blocked by:** 11.

**Status:** ready-for-human

- [x] Terminal viewer (ticket 11 criterion): in a plain terminal REPL on a machine with a GUI (not over SSH), `julia --project=docs`, `include("examples/drawing.jl")`, then `bgplot(motor)`: the drawing should open in the default viewer.
- [x] Plain `connect` to a two-port's `p1` (review of v0.3, item 14): `connect(element.power, tf.p1)` makes a bond pointing into the element, so out of port 1, against "port 1 in". Try it with the examples and decide whether the `TF`/`GY` docstrings need a warning ("reach `p1` through a junction, `Se ⇀ 1 ⇀ TF`"), and whether `PowerPort`'s docstring should say that `e`, `f` are the bond's variables along its half-arrow (MTK's raw connector counts `f` into the component).
- [x] Function naming standard (review of v0.3, item 21): discuss the package's naming conventions, starting with `src/graph.jl`'s internal `source` (renders a `BondGraph` to a format) and `graph` (converts a model, `BondGraph` or plot to a `BondGraph`), proposed as `render` / `to_bondgraph`; extend to the public names if the discussion warrants it.

## Comments

- 2026-10-17: terminal viewer checked by the user: `bgplot(motor)` opens the default viewer.
- 2026-10-18: `p1` check explored in `tmp_dev_ticket12_checks.jl` (gitignored, `julia --project=docs`). Finding: MTK's `connect` is a 0-junction (every bond pointing towards the ports), not a bond; two ports make one bond only if one of them is a sender (negates MTK's flow: `p2`, junction `[-1, :a]`). Elements have only a receiver port, so `connect(se.power, tf.p1)` is `0 ⇀ Se, 0 ⇀ TF`: the Se's flow and power flip silently, and an Sf on `p1` flips the physics (spring q = −2 instead of +2). The sign can't go into `connect` (MTK writes its equations; plain equations instead make MTK zero both flows, q = 0 silently). Prototypes that give `Se ⇀ TF ⇀ I` with correct numbers: E (two-ports take signed members, like junctions) and F (oriented ports: every element has `in` and `out`; a multi-port plain connect then is a correctly signed 0-junction). User sees this as a real limitation (`[-1, x]` with `connect`, submodels); decision open.
- 2026-10-18: `p1` item decided (user): oriented ports, [ADR 0002](../../../docs/adr/0002-oriented-ports.md). Elements get `power` (in) and `pout` (out) plus `e`, `f`; `connect(se.pout, tf.p1)` is `Se ⇀ TF`; junction members `x`, `x => :pout`, `:a`, `:a => :in`; `bondgraph` errors on an element bonded through both ports. Implemented in `.scratch/oriented-ports/` (tickets 01–03). The naming item stays open.
- 2026-10-18: naming decided (user): Julia Base style, lowercase squashed where readable. `bondgraph(model)` and internal `graph(x)` → constructor `BondGraph(x)`; `bgplot` kept; `savebg` → `bgsave`; `bg_equations` → `bgequations`; internal `source` → `render`; `Mass(m, u)` → `Mass(m, v)`; `Spring`/`Damper` default to 1.0. AGENTS.md style line corrected (was "snake_case"). Commit 586ff6e. ADRs and earlier tickets keep the old names as written.
