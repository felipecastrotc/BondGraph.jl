# Choose the dependency and license stance

Type: grilling
Status: deferred
Blocked by: 02

## Question

Should BondGraph depend on full ModelingToolkit (with its AGPL-3.0 pieces) or only on ModelingToolkitBase (MIT, without index reduction/tearing), and does BondGraph's own license (currently Apache-2.0) need to change?

## Deferred (2026-10-14)

The maintainer keeps this an open topic: becoming AGPL in the end is acceptable if required. Until decided:
- **License:** unchanged (currently Apache-2.0; the maintainer said "keep it MIT", to be confirmed which was meant).
- **Dependency (provisional):** depend on full ModelingToolkit (index reduction/tearing needed for higher-index bond graphs; AGPL acceptable). Revisit together with the license before registration / the software paper.
- TODO tracked in `.scratch/roadmap.md`.
