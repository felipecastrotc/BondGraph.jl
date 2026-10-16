# 02: Single-junction mechanical models: 1-DoF and Duffing

**What to build:** A modeller writes the 1-DoF mass–spring–damper and the Duffing oscillator with role-named elements. They can also use the mechanical shorthands `Mass(m, u)`, `Spring(k, x)` and `Damper(c)`, and write a custom nonlinear C with `@oneport`. Both examples reproduce v0.2's golden data.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** ready-for-agent

- [ ] Shorthands: `Mass(m, u)` is `Inertia(I = m, p = m·u)`; `Spring(k, x)` is `Compliance(C = 1/k, q = x)`; `Damper(c)` is `Resistance(R = c)`. v0.2 defaults are kept. They're exported, with docstrings calling them mechanical shorthands.
- [ ] The 1-DoF and Duffing examples are ported with role names. Duffing keeps its linear C and its custom cubic C on one 1-junction. Both pass their golden CSVs.
- [ ] Checks:
  - [ ] `d(p²/2I)/dt = e·f` along a trajectory.
  - [ ] Each shorthand gives the same trajectory as its role-named element.
  - [ ] Two inertias on one 1-junction compile and simulate.
- [ ] Docs examples for the one-ports, the shorthands and `@oneport`.
