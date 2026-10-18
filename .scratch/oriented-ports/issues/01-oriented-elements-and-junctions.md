# 01: Oriented elements and junction members

**What to build:** `@oneport` adds `power`, `pout`, `e`, `f`; the built-in one-ports and adapters write their laws in `e`, `f`. Junctions take `x`, `x => :pout`, `:a`, `:a => :in`; `[-1, x]` errors with the new forms in the message. Convert every test (and regression test) to the new syntax.

**Blocked by:** None.

**Status:** done (feature/mtk-upgrade, uncommitted), user review pending

## Acceptance criteria

- [x] `Pkg.test()` green, regression data unchanged.
- [x] Tests: `Se ⇀ TF ⇀ I` and `Sf ⇀ TF ⇀ C` with plain `connect(x.pout, tf.p1)`; a plain multi-port `connect` of element ports equals `Junction0`; `[-1, x]` throws.
- [x] Docstrings of `PowerPort`, `@oneport`, the elements, junctions and adapters describe the two ports and `e`, `f`.

## Comments

- 2026-10-18 (agent): `bond_ports()` (internal) gives every element `power`, `pout`, `e`, `f` and the three tie equations; `@oneport` and the adapters use it, built-in laws are in `e`, `f`. Junction core: bond flow = sign · core-port flow for every member kind, so the old element/free-port σ split is gone. Member `x.pout` was impossible (a junction owns its elements, and `x.pout` doesn't lead back to `x`), so the form is `x => :pout` (spec, Decisions). Tests converted with a rewriter (`x`→`x => :pout`, `[-1, x]`→`x`, `:a`→`:a => :in`, `[-1, :a]`→`:a`, `x.power.e/f`→`x.e/f`), reviewed by hand. New tests: Se ⇀ TF ⇀ I and Sf ⇀ TF ⇀ C with plain connects, plain multi-port connect = closed form, removed forms throw, unconnected port carries no flow. `Pkg.test()`: 162/162, regression data unchanged.
