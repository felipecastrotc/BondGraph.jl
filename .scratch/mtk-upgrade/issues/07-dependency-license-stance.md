# Choose the dependency and license stance

Type: grilling
Status: resolved
Blocked by: 02

## Question

Should BondGraph depend on full ModelingToolkit (with its AGPL-3.0 pieces) or only on ModelingToolkitBase (MIT, without index reduction/tearing), and does BondGraph's own license (currently Apache-2.0) need to change?

## Deferred (2026-10-14)

The maintainer keeps this an open topic: becoming AGPL in the end is acceptable if required. Until decided:
- **License:** unchanged (currently Apache-2.0; the maintainer said "keep it MIT", to be confirmed which was meant).
- **Dependency (provisional):** depend on full ModelingToolkit (index reduction/tearing needed for higher-index bond graphs; AGPL acceptable). Revisit together with the license before registration / the software paper.
- TODO tracked in `.scratch/roadmap.md`.

## Answer (2026-10-18)

- **Dependency:** full ModelingToolkit. Index reduction and tearing are needed (derivative causality, the motor-pendulum example); ModelingToolkitBase doesn't have them.
- **License:** MIT from 0.3.1 on (maintainer's choice, after comparing MIT, Apache-2.0, MPL-2.0 and AGPL-3.0). MIT is the Julia/SciML norm, e.g. ModelingToolkitStandardLibrary is MIT on full MTK. 0.3.0 and earlier stay Apache-2.0. MTK's AGPL parts still apply to programs that combine them. Open check for the maintainer: whether the university's innovation office (NIT) or the funding agency needs to approve the license.
