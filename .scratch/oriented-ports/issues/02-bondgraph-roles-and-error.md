# 02: `bondgraph` reads port roles and rejects an element bonded through both ports

**What to build:** `bondgraph` treats `pout` as an outgoing port and `power` as incoming, so `connect(se.pout, tf.p1)` draws `Se ⇀ TF`. It raises an error naming the element and both of its bonds when an element uses both ports.

**Blocked by:** 01.

**Status:** done (feature/mtk-upgrade, uncommitted), user review pending

## Acceptance criteria

- [x] `bondgraph` of `connect(se.pout, tf.p1)` + `connect(tf.p2, mass.power)` is `Se ⇀ TF ⇀ I`, no 0 node.
- [x] An element connected through `power` and `pout` throws, and the message names it and both bonds.
- [x] `test_graph.jl` covers both.

## Comments

- 2026-10-18 (agent): element ports registered with roles (`power` :into, `pout` :outof); uses of each element's ports collected from junction members and plain connects; an element with both throws `ArgumentError` naming the element and each bond (e.g. `connect(j₊out, C₊power)`). Tests in `test_graph.jl`.
