# Prototype the bond connection mechanism on MTK 11

Type: prototype
Status: open
Blocked by: 03

## Question

Does the connection mechanism that ticket 03 recommends (native MTK connections, or the BondGraph-owned `expand_bonds` pre-pass from the research note as fallback) actually work on MTK 11? A throwaway prototype should build the RC circuit and the DC motor (gyrator) and reproduce their reference trajectories from `test/regression/reference/`, so the maintainer can react to how model-building code looks with it.
