# 04: Two-ports: `@twoport`, `TF`, `GY` and the DC motor

**What to build:** A modeller uses `TF(m)` and `GY(r)` (port 1 in, port 2 out; Karnopp Eqs. 3.9, 3.11), with a constant or an expression as the modulus. They can write their own two-port with `@twoport` in the bond variables `e1, f1, e2, f2`. The DC motor example reproduces v0.2's golden data.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** ready-for-agent

- [ ] `@twoport` adds `p1` (bond in) and `p2` (bond out) and binds `e1 = p1.e`, `f1 = p1.f`, `e2 = p2.e`, `f2 = −p2.f` before the user's laws.
- [ ] `TF(m)`: `e1 = m e2`, `f2 = m f1`. `GY(r)`: `e1 = r f2`, `e2 = r f1`. Both are written with `@twoport`, and the docstrings note that an expression modulus needs `@named`.
- [ ] The DC motor example is ported and passes its golden CSV.
- [ ] Check: `e1 f1 = e2 f2` along a trajectory for constant and modulated TF and GY.
- [ ] Check: a custom `@twoport` TF matches the built-in one.
- [ ] Docs examples for TF/GY (constant and modulated) and a custom `@twoport`.
