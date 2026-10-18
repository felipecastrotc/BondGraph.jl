# 2. Oriented ports: a bond's direction is in its ports

Date: 2026-10-18
Status: accepted; amends [ADR 0001](0001-junctions-on-native-mtk-connections.md) (element ports, junction members). Implemented on `feature/mtk-upgrade` (`.scratch/oriented-ports/`).

## Context

MTK's `connect` is a 0-junction, not a bond: it equates the efforts and sums the flows to zero, each flow counted into the component that owns the port. Read as a bond graph, every bond of a `connect` points away from it, towards the ports. Two ports therefore make one bond only if one of them negates MTK's flow: a two-port's `p2` (`f2 = −p2.f`) or a junction's outgoing free port. Two ports that both count the flow into themselves make a sign-flipping 2-port 0-junction (Karnopp §3.3, p. 63).

Elements had a single port, `power`, whose bond points into the element. An element could only be at the head of a plain-connect bond, so `connect(se.power, tf.p1)` is `0 ⇀ Se, 0 ⇀ TF`, not `Se ⇀ TF` (ticket 12). It runs without any warning, but the source's flow and power come out negated. With an Sf it changes the physics: a spring behind the TF ends at q = −2 instead of +2. The direction a junction member's `[-1, x]` gives can't be used with `connect` or submodels.

The sign can't go into `connect` itself. MTK writes the connect equations, and changing them means overriding its internals, which ADR 0001 rules out. Writing a bond as plain equations instead of a `connect` doesn't work either: MTK sets the flow of every unconnected port to zero, and the spring above stays at q = 0, again with no error.

## Decision

- **Every element has two power ports: `power`**, whose bond points into the element, **and `pout`**, whose bond points out of it. This covers one-ports, detectors, adapters and `@oneport` elements. The modeller connects one of the two; MTK zeroes the unused one's flow.
- **A bond is a `connect` from an outgoing port to an incoming one.**
  - Outgoing: an element's `pout`, a two-port's `p2`, a junction's outgoing free port.
  - Incoming: an element's `power`, a two-port's `p1`, a junction's incoming free port.

  So `connect(se.pout, tf.p1)` is `Se ⇀ TF`, and `connect(tf.p2, mass.power)` is `TF ⇀ I`.
- **Every element has the variables `e` and `f`.** They are its bond's effort and flow along the half-arrow, whichever port carries the bond. The laws are written in them, and results are read from them (`mass.e`, `se.f`). A port's raw `f` keeps MTK's sign: `se.pout.f` is minus the bond flow.
- **Junction members:** an element `x` is a bond from the junction into `x`; `x => :pout` is a bond from `x` into the junction. A free port `:a` is a bond leaving the junction; `:a => :in` is a bond entering it. The `[-1, x]` form is removed and raises an error that names the new forms. (`x.pout` can't be a member: a junction owns its elements as subsystems, and the connector `x.pout` doesn't lead back to `x`.)
- **Two-ports are unchanged:** `p1` is the bond in, `p2` the bond out.
- **`bondgraph` raises an error** when an element has both of its ports connected, naming the element and both connections. There is no separate compile step: a model that is never drawn isn't checked.

## Considered options

- **Docstring warning only.** It leaves the limitation in place.
- **Each element gets its own default direction for a plain `connect`** (sources out, R/C/I in). The direction would still be invisible in the `connect`, and a source drawn into the system, like Borutzky's `1 → Se:0` (Fig. 4.8), would need a junction.
- **Two-ports take signed members, like junctions** (`TF(se, [-1, mass])`). It fixes two-ports only, and puts the elements inside the TF's namespace.
- **A builder that inserts a hidden 1-junction for `[-1, x]`.** `connect` can only produce an equation, so this needs a replacement for `System`, and it adds nodes the modeller never drew.
- **`connect` argument order as the direction**, as in v0.2. MTK keeps the order, but its equations are symmetric, so the drawing and the equations would disagree.

Prototypes of the two-port-member and oriented-port options gave `Se ⇀ TF ⇀ I` and `Sf ⇀ TF ⇀ C` with the expected numbers (gitignored scratch file `tmp_dev_ticket12_checks.jl`).

## Consequences

- **Breaking change**, made before v0.3 is released:
  - `[-1, x]` is gone.
  - A plain junction member `x` now means a bond into `x`; it used to mean a bond out of `x`. Only this repository uses the v0.3 API.
  - All examples and tests are converted, and the v0.2 regression data must match without being regenerated.
- **A plain `connect` of several element ports is a correctly signed 0-junction**, with no `Junction0` needed (`connect(sf.pout, c.power, r.power)`).
- **The junction core no longer converts signs for element members:** the port the modeller chose carries the direction.
- **The design relies on MTK zeroing an unconnected port's flow.** MTK gives every flow variable at least a one-port connection set (`ModelingToolkitBase`, `src/systems/connectors.jl`, `generate_connection_set!`). This is pinned by a test in `test/test_core.jl`.
- **The port names `power` and `pout`** keep `x.power.e` readable, as "the power port". `pout` was chosen over `power_out` for brevity.
