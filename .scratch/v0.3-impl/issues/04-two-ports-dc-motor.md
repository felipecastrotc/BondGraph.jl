# 04: Two-ports: `@twoport`, `TF`, `GY` and the DC motor

**What to build:** A modeller uses `TF(m)` and `GY(r)` (port 1 in, port 2 out; Karnopp Eqs. 3.9, 3.11), with a constant or an expression as the modulus. They can write their own two-port with `@twoport` in the bond variables `e1, f1, e2, f2`. The DC motor example reproduces v0.2's golden data.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** done (feature/mtk-upgrade)

- [x] `@twoport` adds `p1` (bond in) and `p2` (bond out) and binds `e1 = p1.e`, `f1 = p1.f`, `e2 = p2.e`, `f2 = −p2.f` before the user's laws.
- [x] `TF(m)`: `e1 = m e2`, `f2 = m f1`. `GY(r)`: `e1 = r f2`, `e2 = r f1`. Both are written with `@twoport`, and the docstrings note that an expression modulus needs `@named`.
- [x] The DC motor example is ported and passes its golden CSV.
- [x] Check: `e1 f1 = e2 f2` along a trajectory for constant and modulated TF and GY.
- [x] Check: a custom `@twoport` TF matches the built-in one.
- [x] Docs examples for TF/GY (constant and modulated) and a custom `@twoport`.

## Comments

- 2026-10-15: done. `TF(m)`/`GY(r)`: a non-symbolic modulus becomes a parameter, a symbolic one (`Num`) is used as the modulated ratio. A plain `connect` between a two-port port and an element is a bond pointing into the element (MTK counts flow into each component), so sources reach port 1 through a junction, as in the books. Docs example: `examples/two_ports.jl`; DC motor example ported.
