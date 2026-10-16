# 03: Multi-junction models: 2-DoF and N-DoF

**What to build:** A modeller bonds junctions together through named free ports and a plain `connect`, which acts as a 0-junction whose directions come from the signs at the far ends. The 2-DoF and N-DoF examples reproduce v0.2's golden data. N-DoF is drawn in the textbook orientation `0 → sd`.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 01.

**Status:** ready-for-agent

- [ ] The 2-DoF and N-DoF examples are ported. N-DoF uses `0 → sd` with the coupling springs' initial `q` negated, as in the prototype.
- [ ] The regression helpers map the sign of the N-DoF coupling springs' `q` against the v0.2 reference data. Both examples pass their golden CSVs.
- [ ] Check: a junction-to-junction model built with free ports gives the textbook signs, with power balance at every junction.
- [ ] A docs example on free ports and junction-to-junction bonds.
