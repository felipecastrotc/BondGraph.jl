# 12: Manual checks before the release

**What to build:** Nothing to code. The user checks these by hand, then decides whether anything changes.

**Blocked by:** 11.

**Status:** ready-for-human

- [x] Terminal viewer (ticket 11 criterion): in a plain terminal REPL on a machine with a GUI (not over SSH), `julia --project=docs`, `include("examples/drawing.jl")`, then `bgplot(motor)`: the drawing should open in the default viewer.
- [ ] Plain `connect` to a two-port's `p1` (review of v0.3, item 14): `connect(element.power, tf.p1)` makes a bond pointing into the element, so out of port 1, against "port 1 in". Try it with the examples and decide whether the `TF`/`GY` docstrings need a warning ("reach `p1` through a junction, `Se ⇀ 1 ⇀ TF`"), and whether `PowerPort`'s docstring should say that `e`, `f` are the bond's variables along its half-arrow (MTK's raw connector counts `f` into the component).
- [ ] Function naming standard (review of v0.3, item 21): discuss the package's naming conventions, starting with `src/graph.jl`'s internal `source` (renders a `BondGraph` to a format) and `graph` (converts a model, `BondGraph` or plot to a `BondGraph`), proposed as `render` / `to_bondgraph`; extend to the public names if the discussion warrants it.

## Comments

- 2026-10-17: terminal viewer checked by the user: `bgplot(motor)` opens the default viewer.
