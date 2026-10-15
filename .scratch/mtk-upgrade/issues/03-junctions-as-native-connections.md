# Can junctions be expressed with MTK's native connections?

Type: research
Status: resolved

## Question

On ModelingToolkit v11, can a 0-junction be expressed with MTK's native acausal `connect` semantics (shared effort as the potential, flows summing to zero as a `Flow` variable), and a 1-junction as its dual (shared flow, efforts summing to zero)? How would bond direction / the sign of each bond be expressed (orientation conventions, sign flips in the port, or junction components with signed ports)? What do existing MTK-based bond-graph or port-Hamiltonian efforts (including BondGraphs.jl and ModelingToolkitStandardLibrary) do? List each candidate design with what it costs in ergonomics and what it breaks.

## Answer

Yes. See [research/03-junctions-as-native-connections.md](../research/03-junctions-as-native-connections.md) (cited; each claim marked ran/reasoned/UNVERIFIED). Re-verified by the parent session: RC and DC motor on MTK 11.45.4 match the golden CSVs (max |Δ| ≤ 7e-12).

- **0-junction = native multi-port `connect`** on a `Power` port with effort as the shared (Equality) variable and flow as `[connect = Flow]`; plain `mtkcompile`, no pre-pass, no MTK internals. RC, DC motor, and 2-DoF (forced/unforced) reproduced the golden data.
- **1-junction as a dual connector works but is fragile** (MTK pairs connector variables by position: mixing port kinds compiles silently wrong). Recommended instead: a small 1-junction **component** with one port per bond carrying a ±1 orientation.
- **Bond direction:** `connect` argument order carries no information in MTK 11; flow is "positive into the component". Element bonds therefore always point into the element; only junction-to-junction orientation is real and must be stored as ±1 on the 1-junction port.
- **Two-ports:** TF/GY become ordinary two-port components; GY reproduced the DC motor; a modulated GY reading another element's variable compiles and solves. The `coneqs`/tuple-`@named` machinery goes away.
- **Recommendation for ticket 04 (hybrid):** one `Power` port type; 0-junction = native `connect`; 1-junction = signed-port component; TF/GY as two-port components; thin `Junction0(...)`/`Junction1(...; out = ±1)` wrappers to keep today's feel. Costs: junction-to-junction direction moves to `out = ±1`; generated `a▶b` names and state names change; `generate_graph` must be rebuilt from `connect` equations. N-DoF and Duffing not yet tried.
