# Research 03: Can junctions be expressed with MTK's native connections?

Ticket: `.scratch/mtk-upgrade/issues/03-junctions-as-native-connections.md` · Date: 2026-10-14

## Summary

1. **Yes, on MTK v11.45.4, with no pre-pass and no MTK internals.** A 0-junction is exactly what a native `connect` does to a port whose effort is the potential (`Equality`) and whose flow is `Flow`: efforts are made equal and flows sum to zero. **[ran]**
2. A 1-junction can be the dual `connect` (flow as the potential, effort as `Flow`). That reproduces the DC motor **[ran]**, but every element then needs a port "kind", sources change sign on dual ports, and each 0–1 bond needs an oriented adapter. Mixing the two port kinds is **silently wrong** when the variable names match: MTK pairs connector variables by position. **[ran]**
3. Bond direction **cannot** come from `connect` argument order: MTK v11 gives the same equations for `connect(a,b)` and `connect(b,a)` **[ran]**. The only sign MTK applies is "flow positive into the component" (outer connectors are negated). So element bonds always point into the element, and the orientation of each junction–junction bond has to be explicit data, such as a sign on a 1-junction port.
4. **Recommended for the prototype (ticket 04): the hybrid.** Use one `Power` port type (`e` Equality, `f` Flow). A 0-junction is a native multi-member `connect`. A 1-junction is a small component with one signed port per bond. TF/GY are two-port components on the same `Power` port. With thin `Junction0`/`Junction1` sugar this reproduced the RC, DC motor and 2-DoF (forced and unforced) golden data, max deviation ≤ 3e-11, using plain `mtkcompile` **[ran]**.
5. BondGraphs.jl does almost the same thing. Every bond is a two-member native `connect` of `(E, F[connect=Flow])` ports, and both junctions are explicit components (0: Σ F = 0 with no signs; 1: weights ±1 from bond direction). It has no gyrator. SystemBondGraphs.jl and MTKSL do not use `connect` for junctions; MTKSL has no bond-graph code at all.
6. What changes for users: each bond into a 1-junction needs its own port, there is no "hub" port, and junction–junction orientation is set at construction (`out = ±1`), not by `connect` order. The `coneqs` and tuple-`@named` machinery for mGY/mTF goes away, and state names may change (alias selection).

---

## 0. How it works today (BondGraph.jl v0.2 on MTK 8.71)

- `Power` has `e`, `f`, both tagged `[connect = bg]` plus metadata `(op|j0|j1|tpgy|tptf, bgeffort|bgflow)` (`src/bg.jl:27-36`, `src/types.jl:7-20`).
- `Junction0/1(ps...)` creates one **hub** port `power` and, for each element, `connect(s.power, power)` (sign +1: bond element→junction) or `connect(power, s.power)` (sign −1: junction→element) (`src/bg.jl:159-231`, `src/utils/bg_utils.jl:96-136`). Junction–junction bonds are `connect(j_a.power, j_b.power)` at model level (`examples/mechanical_msd_2dof.jl`).
- The overridden `generate_connection_set` keeps each `connect` as its own **ordered** set (`src/mtk.jl:42-54`). `adjmtx2eqs` turns the order into signs: −1 for "leaving", +1 for "arriving". It applies signs **only to the sums** and forces `s = 1` in the equalities (`src/utils/graph.jl:226-251`). Junction–junction bonds get generated variables `a▶b₊e/f` (`gen_con_name`, `graph.jl:300-342`).
- mGY/mTF have two ports `pin`/`pout` with `e₂ ~ g f₁, e₁ ~ g f₂` and `f₁ r ~ f₂, e₂ r ~ e₁`. Their bonds are emitted as extra `connect`s through `coneqs` or tuple `@named` (`src/bg.jl:254-357`, `src/utils/connection_utils.jl:89-112`).
- All of this depends on method overwrites that are impossible on Julia ≥ 1.10 and on connection internals that were rewritten in MTK v10.4.0 (`.scratch/mtk-upgrade/research-mtk-internals.md` §0, §3).

## 1. What MTK v11 native connections actually do (primary source + runs)

Source: ModelingToolkit v11.45.4, `lib/ModelingToolkitBase/src` (abbreviated **MTKB**, base URL `https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/`). Runs: Julia 1.10.12, ModelingToolkit 11.45.4 / ModelingToolkitBase 1.77.4, ModelingToolkitStandardLibrary 2.29.8, OrdinaryDiffEqTsit5 2.1.5, in a throwaway environment.

| # | Fact | Evidence |
|---|---|---|
| F1 | Only `Equality` (the default), `Flow` and `Stream` (plus causal input/output) are documented connection semantics. `Equality` vars are made equal; `Flow` vars "must sum to zero … assuming appropriate signs are used for current flowing in and out". | MTKB `variables.jl#L141-L170`; docs `docs/src/API/model_building.md#L82-L117` |
| F2 | The sign in a Flow sum depends **only** on inner/outer: `push!(add_buffer, cvert.isouter ? -v : v)`. A connector is "outer" if it is a direct subsystem of the system that contains the `connect`. This is Modelica's "positive into the component" convention, which the docs cite ("similar to … the Modelica specification"). | MTKB `systems/connectors.jl#L839-L852`, `#L334-L368`; `model_building.md#L84` |
| F3 | `connect(a, b, …)` stores its systems in order (`Connection(syss)`), but expansion merges them into connected components of a hypergraph, so the order has **no effect** on the equations. **[ran]** `connect(a,b)` and `connect(b,a)` both give `0 ~ a.f + b.f`. | MTKB `systems/connectors.jl#L25-L39`; `systems/connectiongraph.jl#L392`; probe C |
| F4 | A connector that appears in two `connect` statements ends up in **one** merged set. **[ran]** `connect(a,b), connect(a,x)` → `a.e ~ b.e, a.e ~ x.e, 0 ~ a.f + b.f + x.f`. Any connector with Power semantics is therefore an implicit 0-junction "node". | connectiongraph.jl#L392; probe "merge" |
| F5 | An unconnected inner `Flow` variable gets a singleton set, i.e. `0 ~ f`. **[ran]** An unconnected exposed port therefore carries no power. | MTKB `systems/connectors.jl#L813-L825`; probe D (`0 ~ pout₊f(t)`) |
| F6 | Connected connectors are matched **by position** of their unknowns. The only check is that the *names* are a subset (`any(!in(sys1_dvs_set), …)`) and that the types per position are equal. **[ran]** connecting `Power(e Eq, f Flow)` with `PowerDual(f Eq, e Flow)` compiles **without error** and generates `a₊e ~ b₊f` plus a flow sum mixing `a.f` and `b.e`. When the order is the same and the types are swapped, it errors with "Different types of connectors". | MTKB `systems/connectors.jl#L636-L676`, `#L370-L372`; probe A |
| F7 | A custom `connect = X` type falls into the `else # Equality` branch. MTK has no extension point for junction semantics. | MTKB `systems/connectors.jl#L941-L978`; research-mtk-internals.md §3 |
| F8 | `@connector` warns when #Flow ≠ #regular variables, at **construction time**. A port with two Equality vars (`e`, `f`, no Flow) would warn on every instance. `Power(e, f[Flow])` (1 Flow, 1 regular) is silent. | MTKB `systems/connectors.jl#L148-L176`; `systems/abstractsystem.jl#L3061` |
| F9 | `mtkcompile` works directly on these models. No `expand_bonds` pre-pass is needed because no `connect = bg` variables exist. **[ran]** | runs below |
| F10 | `state_priority` metadata lets the compiler prefer certain unknowns as states. | docs `docs/src/API/variables.md#L303-L307`; MTKB `variables.jl#L60` |

### 1.1 Sign algebra (reasoned, then confirmed by runs)

Every port uses one convention: `power.f` is the flow **into the component** that owns the port (F2). Then:

- **Elements.** I: `D(p) ~ e`, C: `D(q) ~ f`, R: `e ~ R f`. These are the textbook laws with the bond pointing into the element, and they are written once with no sign argument. A source that delivers power has a negative port flow. Today's sign on an element bond (`[-1, R]` vs `R`) has no remaining meaning on a 0-junction. On a 1-junction, flipping **all** signs is invariant, which is why the examples that do this (RC `Junction0(R, C)`, 1-DoF `Junction1(m, d, s, [-1, f])`) still match. A *mixed* sign on a passive element (today: `Junction0(R, [-1, C])` → `f_R − f_C = 0`) means an element that generates power. It cannot be expressed natively, and it is unphysical anyway.
- **0-junction.** `connect(p₁,…,pₙ)` ⇒ `e₁ = … = eₙ`, `Σ fᵢ = 0`. No signs are needed, because with flows "into each component", Σ in = Σ out is the same statement as Σ fᵢ = 0.
- **1-junction.** No native node exists for a shared flow on `Power` ports. Model it as a component with ports pₖ (flow into the junction) and bond signs sₖ (+1 = bond into the junction, −1 = out of it): `sₖ·pₖ.f ~ s₁·p₁.f`, `0 ~ Σ sₖ·pₖ.e`. The two-member `connect` then gives the element `f_elem = −pₖ.f = −sₖ·i`, so every **outgoing** bond (to R, C, I) sees the **same, unsigned** junction flow `i`, and the sign appears only in the effort sum. Today's code gets this behaviour by "disabling the sign in the equality" (`graph.jl:246`; `TODO.md` "Expression signs"). The DC-motor and 2-DoF runs reproduce the golden data with this algebra, which is evidence for ticket 06.
- **Junction–junction bonds.** 0–0 and 1–1 bonds merge consistently in either direction. The orientation of a **0–1 bond** is real information: flipping it flips one term in both junction sums. In a tree that only flips the sign of downstream states (and their initial conditions); in a loop it changes the physics (reasoned). Since `connect` order is lost (F3), it must be stored as the 1-junction port sign.

## 2. What existing projects do

| Project | Junctions | Bond direction / signs | Two-ports | Source |
|---|---|---|---|---|
| **BondGraphs.jl** v0.5.1 (jedforrest, last commit 2025-12-19) | Explicit node components. `EqualEffort`: `0 ~ ΣF`, `E₁ ~ Eₖ` (no weights). `EqualFlow`: `0 ~ Σ wₖEₖ`, `w₁F₁ ~ wₖFₖ`. | Each bond is `connect(src_port, dst_port)` on `MTKPort` with `E(t)` (Equality) and `F(t) [connect = Flow]`. Junction weights are set when the bond is added: −1 if the junction is the source, +1 if it is the destination. | TF only: `E₂ = nE₁, F₁ + nF₂ = 0`. No GY in the standard library. Sf writes `0 ~ fs + F` (delivering source has negative port flow). | [systems.jl#L27-L54](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L27-L54), [#L83-L86](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L83-L86), [#L104-L109](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/systems.jl#L104-L109); [AbstractNode.jl#L84-L104](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/structures/AbstractNode.jl#L84-L104); [graphfunctions.jl#L61-L66](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/graphfunctions.jl#L61-L66); [standard.jl#L82-L114](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/src/libraries/standard.jl#L82-L114) |
| BondGraphs.jl compat | `ModelingToolkit = "9, 10, 11"` but `SymbolicUtils = "3"`, while MTK 11.45.4 needs `SymbolicUtils = "4.31"`. It therefore cannot install alongside current MTK 11 (reasoned from compat bounds; **UNVERIFIED** by resolve). It still uses `ODESystem`/`structural_simplify`. | | | [Project.toml](https://github.com/jedforrest/BondGraphs.jl/blob/v0.5.1/Project.toml); [MTK Project.toml#L131](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/Project.toml) |
| **SystemBondGraphs.jl** (TRACER-LULab, last commit 2024-06-28) | Graph walk: for each 1-/0-junction it emits `0 ~ in_sum − out_sum` and chained equalities from in/out neighbours. **No MTK `connect`.** | From graph edge direction. | Own `ODESystem`s. | [ModelGeneration.jl#L74-L125 @2988bf0](https://github.com/TRACER-LULab/SystemBondGraphs.jl/blob/2988bf08df74e925b269efad245425f9289937d7/src/ModelGeneration.jl#L74-L125) |
| **ModelingToolkitStandardLibrary** v2.29.8 | No bond-graph code: `grep -i "bond\|gyrator\|hamiltonian"` over `src/` and `docs/` finds nothing. Across/through domains: electrical `Pin(v, i[Flow])`, rotational `Flange(phi, tau[Flow])`, translational `MechanicalPort(v, f[Flow])`. | Two-pin `OnePort` (`v = p.v − n.v`, `i = p.i = −n.i`). In bond-graph terms, *parallel* = 0-junction and *series* = 1-junction. Translational uses **velocity** as the potential, so an MTKSL translational node is a mechanical **1**-junction. | `EMF` (electro-mechanical gyrator): `k*w ~ v`, `flange.tau ~ -k*i`. `IdealGear` (TF): `phi_a ~ ratio*phi_b`, `0 ~ ratio*flange_a.tau + flange_b.tau`. The DC-motor tutorial wires `EMF` with `connect`. | [Electrical/utils.jl#L10-L54](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/utils.jl#L10-L54); [ideal_components.jl#L316-L344](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/Analog/ideal_components.jl#L316-L344); [Rotational/components.jl#L245-L270](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Rotational/components.jl#L245-L270); [Translational/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Translational/utils.jl#L10-L16), [components.jl#L76-L93](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Translational/components.jl#L76-L93); [dc_motor_pi.md#L46-L67](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/docs/src/tutorials/dc_motor_pi.md#L46-L67) |
| MTK docs (v11.45.4) | No bond-graph or port-Hamiltonian page: `grep -ril bond docs/ lib/ModelingToolkitBase/src src` finds nothing. | | | MTK repo at tag v11.45.4 |
| Port-Hamiltonian Julia packages (`PortHamiltonianSystems.jl`, `PortHamiltonian.jl`, `PortHamiltonianBenchmarkSystems.jl`) | Found via `gh search repos`. None was inspected further; none appears to build on MTK `connect` (**UNVERIFIED**). | | | GitHub search, 2026-10-14 |

## 3. Candidate designs

Shared prototype code (all of it ran on MTK 11.45.4):

```julia
using ModelingToolkit; using ModelingToolkit: t_nounits as t, D_nounits as D
@connector function Power(; name)            # effort = potential, flow = Flow (into component)
    vars = @variables begin e(t), [guess = 0.0]; f(t), [connect = Flow, guess = 0.0] end
    System(Equation[], t, vars, []; name)
end
I(; name, m)  = (@named power = Power(); @parameters I = m;  System([D(power.f) ~ power.e/I], t, [], [I]; name, systems=[power]))
C(; name, k)  = (@named power = Power(); @parameters C = 1/k; @variables q(t)
                 System([power.e ~ q/C, D(q) ~ power.f], t, [q], [C]; name, systems=[power]))
R(; name, c)  = (@named power = Power(); @parameters R = c;  System([power.e ~ R*power.f], t, [], [R]; name, systems=[power]))
Se(E; name)   = (@named power = Power(); System([power.e ~ E], t, [], []; name, systems=[power]))
function Junction1(signs; name)              # sign +1: bond into junction, -1: out of it
    ps = [Power(name = Symbol(:p, k)) for k in eachindex(signs)]
    System([[signs[k]*ps[k].f ~ signs[1]*ps[1].f for k in 2:length(ps)]; 0 ~ sum(s*p.e for (s, p) in zip(signs, ps))],
           t, [], []; name, systems = ps)
end
GY(; name, g) = (@named p1 = Power(); @named p2 = Power()     # flow into GY at both ports
                 System([p1.e ~ -g*p2.f, p2.e ~ g*p1.f], t, [], []; name, systems=[p1, p2]))
TF(; name, r) = (@named p1 = Power(); @named p2 = Power()     # p1.e*p1.f + p2.e*p2.f == 0 (checked symbolically)
                 System([p1.e ~ r*p2.e, p2.f ~ -r*p1.f], t, [], []; name, systems=[p1, p2]))
```

### D-H. Hybrid: 0-junction = native `connect`, 1-junction = signed component **(recommended)**

RC **[ran, Δ = 1.1e-16]**:
```julia
@named Cc = C(k = 1.0); @named Rr = R(c = 1.0)     # (Cc/Rr only to avoid shadowing the constructors)
@named model = System([connect(Rr.power, Cc.power)], t; systems = [Rr, Cc])   # or sugar: Junction0(Rr, Cc)
sys = mtkcompile(model)                                                        # D(Cc₊q) ~ Cc₊power₊f
```
DC motor **[ran, Δ = 7.2e-13]** (no sugar shown; the sugar version below is design D-E):
```julia
@named je = Junction1([+1, -1, -1, -1]); @named jm = Junction1([+1, -1, -1, +1]); @named gy = GY(g = 0.01)
eqs = [connect(V.power, je.p1), connect(R.power, je.p2), connect(L.power, je.p3), connect(je.p4, gy.p1),
       connect(τ.power, jm.p1), connect(c.power, jm.p2), connect(J.power, jm.p3), connect(gy.p2, jm.p4)]
sys = mtkcompile(System(eqs, t; systems = [V, R, L, τ, c, J, je, jm, gy], name = :mdl))
```
2-DoF with today-style nested sugar **[ran, unforced Δ = 3.3e-16, forced (`Se(sin(t))`) Δ = 2.9e-11]**. The sugar wraps a 1-junction and exposes one outer `power` port. It plays the role of today's hub, but its orientation is fixed by `out`:
```julia
@named m1_j1 = Junction1([-1, m], [-1, d], [-1, s]; out = -1)   # bond m1_j1 → 0
@named m2_j1 = Junction1([-1, m], f; out = +1)                  # bond 0 → m2_j1
@named m2_sd = Junction1([-1, s], [-1, d]; out = +1)            # bond 0 → m2_sd
@named mdl = System([connect(m1_j1.power, m2_j1.power, m2_sd.power)], t;   # the 0-junction
                    systems = [m1_j1, m2_j1, m2_sd])
```
- **Ergonomics vs today.** `Junction0(R, C)` and `Junction1(V, [-1, R], [-1, L])` survive as thin wrappers, built in about 10 lines of sugar. Users can also write raw `connect`, and any MTK user can read the result. A bond between junctions is a plain `connect` of their `power` ports. Its orientation moves from `connect` argument order to `out = ±1` on the 1-junction. mGY/mTF become ordinary components joined by `connect`, so `coneqs`, tuple `@named` and `gen_tp_con!` disappear. **Modulated** ratios work: `@named gy = GY(g = 0.01 + 0.001*J.power.f)` (sibling variable through `@named`'s `ParentScope`, MTK NEWS L478-484) compiled and solved **[ran]**.
- **What breaks.** The hub port semantics: a 1-junction needs one port per bond, and a raw `connect(j1.power, a, b)` on a 1-junction port builds a 0-junction. Signs on 0-junction arguments become meaningless (warn or reject them). Generated `a▶b₊e/f` names are gone: a bond's variables are the port variables. `generate_graph` must be rebuilt from `connect` equations; a `connect` with more than 2 members is a 0-node. State names change because tearing picks aliases (`jm₊c₊power₊f` instead of `jm₊J₊power₊f` **[ran]**). Fix this with an explicit momentum state plus `state_priority` (F10), or remap names in `test/regression/helpers.jl`. `Sf` must write `power.f ~ -expr` to keep "delivers flow" (as BondGraphs.jl does).
- **MTK internals:** none. Public API only: `@connector`, `Flow`, `connect`, `System`, `mtkcompile`. `mtkcompile` on a raw model is correct.

### D-E. Explicit components for both junctions (BondGraphs.jl style)

Same as D-H, but the 0-junction is also a component (`Junction0(n)`: `pₖ.e ~ p₁.e`, `0 ~ Σ pₖ.f`), and every bond is a two-member `connect`. RC **[ran, Δ = 1.1e-16]**: `connect(R.power, j0.p1), connect(C.power, j0.p2)`. DC motor with today-style sugar **[ran, Δ = 6.2e-12]**:
```julia
@named je = Junction1(V, [-1, R], [-1, L], -1)    # trailing Int = free port p4, bond je → gy
@named jm = Junction1(τ, [-1, c], [-1, J], +1)    # free port p4, bond gy → jm
@named gy = GY(g = 0.01)
@named mdl = System([connect(je.j.p4, gy.p1), connect(gy.p2, jm.j.p4)], t; systems = [je, gy, jm])
```
- **Ergonomics:** the same as D-H. Every junction is a named object, which makes graph drawing and "which junction is this" introspection trivial. It has a few more alias equations, which `mtkcompile` removes.
- **What breaks / internals:** the same as D-H. One extra guard is needed: a user who writes a three-member `connect` on these ports silently creates an extra 0-junction (F4).
- Variant D-E′ (research-mtk-internals §3 option C): ports with `e`, `f` both Equality, so a bond is an identification and every sign lives in the junction. This triggers the F8 warning on every port, and a three-member `connect` makes `e` *and* `f` all equal, which is meaningless. It is worse than D-E with no benefit (reasoned).

### D-N. Fully native: 0 = `connect` on `Power`, 1 = `connect` on a dual port

`PowerDual`: `f` Equality, `e` Flow. DC motor **[ran, Δ = 7.2e-13]**:
```julia
@named L = I(m = 0.5, kind = :dual); @named R = R(c = 1.0, kind = :dual); @named V = Se(12.0, kind = :dual) # Se writes e ~ -E
...
@named gy = GYdual(g = 0.01)                       # p1.e ~ g*p2.f, p2.e ~ -g*p1.f
eqs = [connect(V.power, R.power, L.power, gy.p1),  # electrical 1-junction
       connect(τ.power, c.power, J.power, gy.p2)]  # mechanical 1-junction
```
- **Ergonomics:** the most "MTK-pure" equations for pure 0-only or 1-only models. But every element constructor needs a `kind`, chosen *before* the element is placed on a junction. Sources change sign with the kind, because on a dual port the effort is the oriented, summed variable. Every 0–1 bond needs an oriented adapter component (e.g. `f_a ~ ±f_b, e_b ~ ∓e_a`; reasoned, **not run**). Alternating 0/1 graphs (the mechanical N-DoF examples) need one adapter per junction–junction bond.
- **What breaks:** the same as D-H, plus the element API. **Hazard (F6, ran):** `connect(Power, PowerDual)` with variables named `e`/`f` in both is accepted silently and wrong. Fix it by giving the dual port different variable names, which turns the mistake into a `connection_error`.
- **MTK internals:** none.

### D-P. Two-pin across/through (MTKSL style)

Elements become `OnePort`s with `p`/`n` pins. A 0-junction is a parallel connection and a 1-junction is a series connection; a GY is an `EMF`. RC: `connect(R.p, C.p); connect(R.n, C.n, ground.g)`. DC motor: the MTKSL tutorial wiring (`dc_motor_pi.md#L46-L67`). **Not run.**
- This is a circuit/Modelica API, not a bond graph. It needs a ground reference, the series/parallel duality depends on the domain (translational velocity is the potential, so a mechanical 1-junction becomes a node), and junction structure is no longer explicit. It breaks the dialect's purpose and does not fit the invariant in `map.md`.
- It is still useful as the **coupling target**. Since `Power` has variables `e, f` and `Pin` has `v, i`, a direct `connect` is rejected by the name check (F6, reasoned from source). Coupling needs adapter components (an MTKSL translational node needs a dual or unit-GY adapter). That belongs to the map's "coupling" item.

### D-X. Keep today's API with a BondGraph pre-pass (`expand_bonds`)

Fallback from `research-mtk-internals.md` §3 approach A. It keeps `connect` order as bond direction by reading `Connection.systems` before MTK merges the sets (F3). Costs: `mtkcompile(raw_model)` silently produces Equality physics (F7). It relies on argument order inside `Connection`, which is stored but not promised by the docs. It rewrites `graph.jl`/`connection_utils.jl`. **Not needed** now that D-H/D-E work natively.

### D-C. Custom connection type (`connect = bg`)

Not possible without MTK internals: F7 plus research-mtk-internals.md §3. Rejected.

## 4. Two-ports (TF, GY, mTF, mGY)

- All designs except D-N use two `Power` ports with flow into the two-port at both ports. Power conservation is then `p1.e·p1.f + p2.e·p2.f = 0`:
  GY `p1.e ~ -g·p2.f, p2.e ~ g·p1.f` (matches the DC-motor golden data, **ran**); TF `p1.e ~ r·p2.e, p2.f ~ -r·p1.f` (power balance checked symbolically, **ran**; no golden TF example). These have the same structure as BondGraphs.jl `TF` and MTKSL `EMF`/`IdealGear` (§2).
- Today's `e₁ ~ g f₂` uses the "power in at 1, out at 2" convention. With into-component ports, one sign moves into the constitutive law. This is a convention change for ticket 06, not a physics change.
- Modulated forms only change the expression for `g`/`r`. A ratio that reads a sibling's variable works through `@named`'s `ParentScope` **[ran]**. A causal `RealInput` modulation port (MTKSL Blocks) would be the more MTK-idiomatic alternative (**UNVERIFIED**, not run).
- Two-ports sit on 1-junctions through the signed port (D-H/D-E) or directly in a 0-`connect`. In D-N they need dual ports.

## 5. Can one element port take part in both a 0- and a 1-junction?

- In a bond graph each element port carries **one** bond, so it attaches to one junction. What alternates is the junction type, linked by junction–junction bonds.
- Natively, a port belongs to exactly **one** merged connection set (F4), and a set has one shared variable. A `Power` port therefore cannot join a "shared-flow" set.
- In D-H/D-E there is one port type, and the 1-junction is a component whose ports join ordinary sets. A two-member set is a bond, and a larger one is a 0-junction. 0–1–0 alternation works (2-DoF **ran**).
- In D-N the alternation needs adapters, and connecting the two port kinds without one is silently wrong unless the variable names differ (F6 **ran**).

## 6. Comparison

| | D-H hybrid | D-E explicit (BondGraphs.jl) | D-N dual native | D-P two-pin | D-X pre-pass |
|---|---|---|---|---|---|
| RC / DC motor / 2-DoF reproduce golden | ✔ / ✔ / ✔ (ran) | ✔ / ✔ / – (ran) | ✔ / ✔ / – (ran) | not run | not run (old algorithm) |
| `mtkcompile(raw)` correct | ✔ | ✔ | ✔ | ✔ | ✘ silently wrong |
| MTK internals | none | none | none | none | `Connection.systems` order |
| Element API change | none (single port type) | none | `kind` on every element; sources flip | rewrite to pins | none |
| Bond direction expressed by | J1 port sign / `out` | J1 port sign | element kind + adapter orientation | pin polarity | `connect` order (today) |
| Main user-visible loss | `connect` order no longer means direction; no 1-junction hub | same + more alias eqs | element ergonomics, adapters | not a bond graph | footgun with plain MTK |

## 7. Recommendation for ticket 04 (prototype)

Try **D-H** first, keeping D-E's explicit `Junction0` component as a one-line switch if graph introspection wants named 0-nodes. Concretely:
1. One connector: `Power` with `e` (Equality) and `f [connect = Flow]`, plus guesses only (no initial values). Flow means power into the component.
2. One-ports written once in the "into element" convention. `Sf` writes `f ~ -expr`. Add momentum and displacement states with `state_priority` so states stay `I.p` / `C.q`.
3. `Junction0(args...)` = wrapper that emits one multi-member `connect` and exposes `power`. `Junction1(args...; out)` = signed component plus wrapper exposing `power`. Accept and warn about signs on 0-junction arguments.
4. `GY`/`TF`/`mGY`/`mTF` as two-`Power`-port components connected with `connect`. Drop `coneqs` and tuple `@named`.
5. Port all six examples and compare with `test/regression/reference/*.csv`, remapping state names in `helpers.jl` if needed. N-DoF and Duffing were **not** tried here, and the N-DoF "merged non-adjacent elements" oddity (`map.md`) may surface as a real difference.

Open for ticket 06 (sign convention): the runs support "signs only in the sums when element-side variables are read, into-element port convention, explicit orientation only on 0–1 bonds". Whether to expose `out = ±1` or a `bond(a, b)` helper that infers it is an API decision.

## Method and provenance

- **Read:** MTK v11.45.4 source (`lib/ModelingToolkitBase/src/systems/connectors.jl`, `connectiongraph.jl`, `variables.jl`) and docs (`docs/src/API/model_building.md`, `variables.md`, tutorials); ModelingToolkitStandardLibrary v2.29.8 (fresh clone at tag); BondGraphs.jl v0.5.1 (`eeac405`, latest commit on master per `gh api`); SystemBondGraphs.jl `2988bf0`; this repo's `src/`, `examples/`, `GLOSSARY.md`, `.scratch/mtk-upgrade/*`.
- **Ran:** throwaway scripts in the session scratchpad (`…/scratchpad/mtk11-env/junctions.jl`, `probeA.jl`, `probeB.jl`, `probe2dof.jl`), Julia 1.10.12, separate environment (MTK 11.45.4, MTKBase 1.77.4, MTKSL 2.29.8, OrdinaryDiffEqTsit5 2.1.5). Δ = max |simulated − golden| over the 21 sample times, Tsit5 with abstol = reltol = 1e-10, same as `test/regression/helpers.jl`. The scratchpad is session-local; the code needed to reproduce the runs is in §3.
- **Reasoned, not run:** D-N adapters for 0–1 bonds; D-P; physical meaning of flipping 0–1 bonds in loops; `Power`↔`Pin` rejection; BondGraphs.jl incompatibility with MTK 11 (compat bounds only); `RealInput` modulation.
