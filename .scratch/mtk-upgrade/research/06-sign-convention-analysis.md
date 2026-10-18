# Sign convention and bond direction: evidence for ticket 06

## Summary
1. v0.2 encodes direction in `connect` argument order; junction **sums** use incoming `+` and outgoing `−`, while shared-variable **equalities are unsigned**.
2. This is the rule stated by [20-sim](https://www.20sim.com/webhelp/modeling_tutorial_bond_graphs_frombondgraphtoequations.php); applying the sign to the DC motor's flow equality changes the back EMF from `−gω` to `+gω`.
3. The six regression equation files match the equations of their **code-built** graphs after algebraic cancellation; the N-DoF file contains nonlocal, cancelling spring terms.
4. Two drawings disagree with their code: DC motor torque arrow; 2-DoF force and ground-element placement. N-DoF has no dedicated TikZ graph.
5. The accepted v0.3 hybrid can express all six code-built models with unsigned equalities, explicit 1-junction orientation, and MTK-native 0-junctions; the prototype reports numerical agreement.
6. This is a source-and-algebra review. No Julia execution or regeneration was performed here.

## 1. What the current implementation actually does

**Encoding and path to equations.** `Power` tags both `e` and `f` with custom `connect = bg` and separate `(element type, bgeffort|bgflow)` metadata (`src/bg.jl:27-35`, `src/types.jl:7-20`). A bare junction argument gets sign `+1`; `[number, system]` gets `sign(number)` (`src/utils/bg_utils.jl:96-135`). Thus `[-1, s]` reverses a bond; it does not negate a constitutive parameter or port variable. Zero is accepted by `flatinput` but matches neither connection branch, so that argument silently creates **no bond** (`src/bg.jl:172-177,220-225`). `Junction0` and `Junction1` both make one hub `power`: `+1` emits `connect(element.power, junction.power)`, `−1` emits the reverse (`src/bg.jl:159-182,207-230`). Model-level `connect(a.power,b.power)` has its arrow from its first argument to its second. Neither junction accepts a `couple` keyword today; the docstring mentions it but the signature does not (`src/bg.jl:138-159,185-207`). The old `couple` namespace code is commented out (`src/mtk.jl:95-161`). Metadata identifies which equation rule applies; it carries no arrow bit.

The override of `generate_connection_set` extracts `bg` sets and merges the other sets; `expand_connections` passes the resulting domain sets to the custom graph generator (`src/mtk.jl:25-54`, `src/utils/connection_utils.jl:18-38`). `csets2adjmtx` puts an edge from the first entry of each set to every later entry (`src/utils/graph.jl:144-183`). At node `i`, an outgoing edge gets `−1` and an incoming edge `+1` (`src/utils/graph.jl:222-233`). `gen_con_name` gives a junction-to-junction bond an oriented `a▶b₊e/f` alias; its docstring and branch make positive sign mean `con_var → base_var` (`src/utils/graph.jl:235-242,295-342`). The signed expression is then used in the sum **only**: `s = apply_equalityeqs ? 1 : s` (`src/utils/graph.jl:243-259`). `equalityeqs` equates the first expression to each other, without a sign of its own; `sumvar` creates `0 ~ sum(expressions)` (`src/utils/bg_utils.jl:26-41,67-74`). The local hub variable itself is omitted from both generated rules (`src/utils/graph.jl:235-259`). The implementation's comments explicitly say it disabled the sign for equalities (`src/utils/graph.jl:243-249`); `TODO.md:5-15` records the earlier doubt, not the current branch's behavior.

| Junction | Signed conservation equation at its hub | Unsigned equality | Source |
|---|---|---|---|
| 0 | `0 = Σ(incoming f) − Σ(outgoing f)` | `e₁ = e₂ = …` | `src/utils/graph.jl:209-259` |
| 1 | `0 = Σ(incoming e) − Σ(outgoing e)` | `f₁ = f₂ = …` | `src/utils/graph.jl:209-259` |

These are equations for *bond variables* obtained from the adjacency graph. Element constitutive laws remain `D(f_I)=e_I/I`, `e_C=q/C, D(q)=f_C`, `e_R=R f_R`, and `e_Se=source` (`src/bg.jl:58-63,379-384,406-413,434-440`). Consequently a reversed arrow at a one-port cannot safely be interpreted as simply negating that element's `f` without translating its local power convention too. `mGY` has `e₂=g f₁, e₁=g f₂`; `mTF` has `r f₁=f₂, r e₂=e₁` (`src/bg.jl:254-268,317-331`). `gen_tp_con!` connects the first subsystem *to* `pin`, and `pout` *to* the second subsystem (`src/utils/connection_utils.jl:89-111`). The two-port equations have no junction sign parameter. Their power interpretation depends on those pin/pout arrows; port 1 receives and port 2 sends equal power for these equations. With `mGY(je,jm)`, that is `je → gy → jm` (`examples/electrical_dc_motor.jl:76-90`).

The external reference agrees with the implemented **junction algebra**: [20-sim, “From Bond Graph to Equations,” junction step and DC motor example](https://www.20sim.com/webhelp/modeling_tutorial_bond_graphs_frombondgraphtoequations.php) states shared effort at 0, shared flow at 1, and adds/subtracts only the summed variable according to arrow direction; its motor writes `i_in=i_r=i_l=i_m` and `u_in−u_r−u_l−u_m=0`. Its sentence about efforts at a 1-junction says “towards the 0 junction,” apparently a typo in that 1-junction sentence; its displayed 1-junction equation removes the ambiguity. The TikZ helper draws `\drawbg{A}{B}` with the half-arrow at **B** (`docs/src/misc/tikz/base_tikzset.tex:31-55`). Hence the direction lists below come from the actual source, not the positions of the labels.

## 2. Hand derivations against the saved simplified equations

Notation: `v=ẋ`, `i` is armature current, `ω` rotor speed, `q` is the displacement/charge stored in a C element, `k=1/C`, and each `c`/`R` is the corresponding damping/resistance coefficient. The arrows set the *junction sum* sign; equalities carry no arrow sign. When the local one-port flow in `src/bg.jl` points opposite a drawn bond (notably RC), I explicitly translate it to the flow along that bond. The text files are **simplified v0.2 outputs**, not unsimplified junction equations or a fresh run (`.scratch/golden-tests/report.md:3-8`). A match below means algebraic equality of the stated code-built model and the saved expression, not independent numerical validation.

### Electrical RC

Both arrows are `R → 0` and `C → 0` (`docs/src/misc/tikz/bg_rc.tex:11-20`); code passes both as default `+1` (`examples/electrical_rc.jl:31-38`). Let `f_R^b,f_C^b` be flows **along the drawn arrows**. Standard 0-junction: `e_R=e_C=V_C`, `f_R^b+f_C^b=0`. Since both arrows leave the passive elements, their local into-element flows are `−f_R^b,−f_C^b`: `V_C=R(−f_R^b)=q/C` and `q̇=−f_C^b`. Therefore `q̇=−q/(CR)`, equivalent to the prose `C V̇+V/R=0` (`examples/electrical_rc.jl:7-13`). Saved line 1 is exactly `D(C.q)=−C.q/(C.C*R.R)` (`test/regression/reference/electrical_rc_equations.txt:1`). **Match:** shared voltage, resistor current, capacitor charge derivative, and their negative sign. Caveat: the present implementation's raw `f_R,f_C` are local element variables, so the saved equation alone does not prove that the TikZ arrow is honored as a power sign at each passive port.

### Electrical DC motor

Drawn electrical arrows: `V → je`, `je → L,R,gy`; mechanical: `gy → jm`, `jm → J,c,τ` (`docs/src/misc/tikz/bg_dcmotor.tex:11-43`). The *code* instead has `τ → jm` through the default `+1`, and `jm → J,c` (`examples/electrical_dc_motor.jl:64-82`). For the code-built graph, unsigned flow equalities give `f_V=f_R=f_L=f_{gy,1}=i` and `f_{gy,2}=f_τ=f_c=f_J=ω`. The signed effort sums and gyrator laws are

```text
je: V − R i − L i̇ − e_gy,1 = 0;       e_gy,1 = g ω
jm: τ + e_gy,2 − c ω − J ω̇ = 0;       e_gy,2 = g i
```

Thus `L i̇=V−Ri−gω`, `J ω̇=τ+gi−cω`, exactly the physics prose and direct model (`examples/electrical_dc_motor.jl:13-29,147-153`) and the four source, friction, back EMF, and motor torque terms in saved lines 1–2 (`test/regression/reference/electrical_dc_motor_equations.txt:1-2`). **Match for the code-built graph.** **Mismatch for the drawing:** with its `jm → τ` arrow and the same positive effort `τ`, the standard junction sum yields `J ω̇=gi−cω−τ`. The source polarity or that one arrow must be changed to depict the positive applied torque of the prose/code. This drawing mismatch is separate from whether equality signs are disabled.

### Mechanical 1-DoF

The drawing has `F → 1 → m,d,s` (`docs/src/misc/tikz/bg_1dof.tex:11-27`). Its 1-junction is `v_m=v_d=v_s=v_F=v`, `F−m v̇−cv−q/C=0`, `q̇=v`. The example passes `m,d,s` as `+1` and `f` as `−1`, reversing **all** drawn arrows (`examples/mechanical_msd_1dof.jl:49-61`). Reversing the whole junction sum multiplies it by `−1`, so its physics is unchanged when equalities stay unsigned. Saved lines 1–2 give `v̇=(F−q/C−cv)/m`, `q̇=v` (`test/regression/reference/mechanical_msd_1dof_equations.txt:1-2`). **Match term by term:** force `+`, spring `−`, damper `−`, shared velocity `+`. The drawing and code arrow choices nevertheless do not agree literally.

### Mechanical 2-DoF

For the **code-built** graph, define `v₁=m1_j1.m.f`, `v₂=m2_j1.m.f`, `q₁=m1_j1.s.q`, and `q_b=m2_sd.s.q`. `m1_j1 → 0 → m2_j1,m2_sd`, with all passive one-ports sent *from* their 1-junctions; the source enters `m2_j1` (`examples/mechanical_msd_2dof.jl:66-90`). The unsigned flow equalities and the 0-junction flow sum give `q̇₁=v₁`, `q̇_b=v₁−v₂`. The signed effort sums give

```text
m₁ v̇₁ = −k₁ q₁ − c₁ v₁ − k_b q_b − c_b(v₁−v₂)
m₂ v̇₂ = +F + k_b q_b + c_b(v₁−v₂)
```

Saved line 1 contains `−q₁/C₁`, `−q_b/C_b`, `−c₁v₁`, and `−c_b(v₁−v₂)`; line 2 is `q̇₁=v₁`; line 3 contains `+F`, `+q_b/C_b`, and `+c_b(v₁−v₂)`; line 4 is `q̇_b=v₁−v₂` (`test/regression/reference/mechanical_msd_2dof_equations.txt:1-4`). **All four code-built equations match term by term.**

The **drawing** differs: `F` enters the **first** mass junction, while that junction has only its mass; the **second** junction has the ground spring and damper (`docs/src/misc/tikz/bg_2dof.tex:11-52`). Applying the same standard rules to that drawing, with `q_b` for its intermass branch and `q₂` for its second-mass ground spring, gives `q̇_b=v₁−v₂`, `q̇₂=v₂`, `m₁v̇₁=F−k_bq_b−c_b(v₁−v₂)`, `m₂v̇₂=k_bq_b+c_b(v₁−v₂)−k₂q₂−c₂v₂`. **Mismatch with the saved file:** `+F`, `−k_ground q_ground`, and `−c_ground v` are attached to opposite mass rows, and the ground-spring derivative is `v₂` in the drawing versus `v₁` in the code. The bridge spring/damper terms and `q̇_b` match. The prose's displayed equation is explicitly marked TODO and is not a usable two-mass specification (`examples/mechanical_msd_2dof.jl:14-23`).

### Mechanical N-DoF (`n=4`)

There is **no** `docs/src/misc/tikz/bg_*ndof*.tex`: the example itself marks its bond graph image and governing equation as TODO (`examples/mechanical_msd_ndof.jl:5-20`). The reproducible graph is the construction loop (`examples/mechanical_msd_ndof.jl:46-81`): mass 1 has ground `s,d`; each `i=2…4` has a spring/damper branch `sd_i → 0_i`, a mass branch `0_i → 1_i → m_i`, and the preceding mass junction points into `0_i`. For `v_i=m_i.f`, `q₁` the ground spring displacement, and `q_i` the branch spring displacement, standard unsigned flow equality plus signed 0-junction flow conservation gives `q̇₁=v₁`, `q̇_i=v_i−v_{i−1}` for `i=2,3,4`. Define `B_i=k_iq_i+c_i(v_i−v_{i−1})` and `G₁=k₁q₁+c₁v₁`. The hand-derived chain is

```text
m₁ v̇₁ = −G₁ + B₂       q̇₁ = v₁
m₂ v̇₂ = −B₂ + B₃       q̇₂ = v₂−v₁
m₃ v̇₃ = −B₃ + B₄       q̇₃ = v₃−v₂
m₄ v̇₄ = −B₄            q̇₄ = v₄−v₃
```

Saved lines 2,3,5,7 match the four `q̇` equations exactly. Saved line 1 has `−k₁q₁−c₁v₁ +k₂q₂+c₂(v₂−v₁)` **plus** `−k₃q₃+k₃q₃ +k₄q₄−k₄q₄`; line 4 has `−k₂q₂−c₂(v₂−v₁)+k₃q₃+c₃(v₃−v₂)` **plus** `+k₄q₄−k₄q₄`. Lines 6 and 8 are exactly `−B₃+B₄` and `−B₄`, respectively (`test/regression/reference/mechanical_msd_ndof_equations.txt:1-8`). **Match after exact cancellation, mismatch in unsimplified locality:** downstream springs appear as paired terms where the hand-derived local junction has no such bond. There is no claim here that a nonexistent N-DoF drawing independently confirms the arrows. The prose's displayed N-DoF formula is also marked TODO (`examples/mechanical_msd_ndof.jl:8-13`).

### Mechanical Duffing

The drawing has `F → 1 → m,d,s,s3` (`docs/src/misc/tikz/bg_duffing.tex:12-33`), exactly the signs in `Junction1([-1,m],[-1,d],[-1,s],[-1,s3],f)` (`examples/mechanical_duffing.jl:86-99`). Equality gives a common velocity `v`. The two C elements have independent states `q,q₃`, both with derivative `v`; their efforts are `kq` and `k₃q₃³` (`examples/mechanical_duffing.jl:5-16,55-67`, `src/bg.jl:406-413`). Effort conservation gives `m v̇=F−cv−kq−k₃q₃³`, `q̇=v`, `q̇₃=v`. Saved lines 1–3 have precisely `+F`, `−d.R*v`, `−s.q/s.C`, `−s3.q³/s3.C`, and the two `+v` derivatives (`test/regression/reference/mechanical_duffing_equations.txt:1-3`). **Match term by term.** The two q states are equal only if initialized equally; the regression report confirms that condition for the redundant state (`.scratch/golden-tests/report.md:26`).

## 3. DC motor counterfactual: sign in the equality

Use the *code* directions `V → je → R,L,gy.pin` and `gy.pout,τ → jm → c,J` (`examples/electrical_dc_motor.jl:76-82`, `src/utils/connection_utils.jl:89-110`). Anchor state signs to the saved names `i=L.power.f` and `ω=J.power.f`. Keeping the source efforts `V,τ` positive:

| Flow equality rule | Implied two-port flows | Two-port efforts from current `mGY` | Resulting differential equations |
|---|---|---|---|
| **Current, unsigned:** `f_V=f_R=f_L=f_pin=i`; `f_pout=f_τ=f_c=f_J=ω` | `f_pin=i`, `f_pout=ω` | `e_pin=gω`, `e_pout=gi` | `L i̇=V−Ri−gω`; `J ω̇=τ+gi−cω` |
| **Counterfactual, sign also in equality:** `f_V=−f_R=−f_L=−f_pin`; `f_pout=f_τ=−f_c=−f_J` | `f_pin=i`, `f_pout=−ω` | `e_pin=−gω`, `e_pout=gi` | `L i̇=V−Ri+gω`; `J ω̇=τ+gi−cω` |

The counterfactual changes the **back-EMF term** alone with these fixed state definitions and the present pin/pout law. It disagrees with both the direct equations (`examples/electrical_dc_motor.jl:149-153`) and the saved first equation (`test/regression/reference/electrical_dc_motor_equations.txt:1`). This calculation is about changing `src/utils/graph.jl:249` without redesigning `mGY`; a simultaneous two-port polarity change could compensate algebraically, but would be a different convention and needs a separate power-balance specification. The current row matches the regression report's physical interpretation (`.scratch/golden-tests/report.md:23-25`).

## 4. Why N-DoF produces cancelling nonlocal terms

The nearest-neighbor graph in `examples/mechanical_msd_ndof.jl:58-71` has no physical `m₁—s₃` or `m₂—s₄` bond. Yet the saved m₁ row contains `−q₃/C₃+q₃/C₃+q₄/C₄−q₄/C₄`, and the m₂ row has `+q₄/C₄−q₄/C₄` (`test/regression/reference/mechanical_msd_ndof_equations.txt:1,4`). A **wrong sign on an equality is not needed** to explain this: the current generator explicitly removes equality signs (`src/utils/graph.jl:243-255`). The plausible structural cause is the legacy connection-set-to-adjacency path: the N-DoF loop reuses each 1-junction hub in bonds to adjacent 0-junctions, the override merges domain sets (`src/mtk.jl:42-53`), the graph routine makes a star from the first member of each resulting set (`src/utils/graph.jl:144-183`), and it sums every adjacent signed expression (`src/utils/graph.jl:222-259`). The original 0-junction locality is not kept as an explicit, isolated component equation. A downstream spring can therefore be visited along two oppositely oriented graph incidences; after alias elimination, these show up as `+q/C−q/C`. The regression report independently notes merged non-adjacent 0-junction sets (`.scratch/golden-tests/report.md:23-26`). **Limit:** without inspecting the actual unsimplified N-DoF connection sets, the precise MTK merge step responsible for each pair is unproven. The cancellation in the saved text and the locality violation are certain; the exact internal path is an inference. The accepted MTK 11 prototype reports that the anomaly disappears (`.scratch/mtk-upgrade/issues/04-prototype-connection-mechanism.md:23-25`).

## 5. Translation to the proposed v0.3 hybrid

Research 03 defines one `Power(e: Equality, f: Flow)` with `f` positive **into its owning component**, a native multi-member `connect` for 0-junctions, and a component with signed ports for 1-junctions (`.scratch/mtk-upgrade/research/03-junctions-as-native-connections.md:7-12,41-48,63-84`). A junction port sign `s=+1` means its bond points **into the 1-junction**, `s=−1` points **out**. Its equations are `s_k p_k.f=s_1 p_1.f` and `0=Σs_k p_k.e`; the native two-port `connect` adds `f_element=−p_k.f`, so the element's local flow stays positive into that element and the *junction flow equality is unsigned in the common physical flow* (`research/03-junctions-as-native-connections.md:43-48,76-84`). A 0-junction is `e₁=…=eₙ, Σf_port=0` with no per-port sign. MTK 11 ignores `connect` argument order (`research/03-junctions-as-native-connections.md:28-34`), so junction-to-junction direction moves to the 1-junction port sign. A delivering flow source needs `power.f=−Sf_value`; GY becomes `e₁=−g f₂, e₂=g f₁`, and TF `e₁=r e₂, f₂=−r f₁`, since **both** two-port flows now point into the two-port (`research/03-junctions-as-native-connections.md:43-48,75-84,155-160`). Those laws satisfy `e₁f₁+e₂f₂=0`. GY signs cannot simply be copied from v0.2.

| Example | Translate today's signs and orientation | Hand-equation check under the hybrid |
|---|---|---|
| RC | Drop the two 0-junction `+1` signs; `connect(R.power,C.power)`; both element ports use into-element flow. The drawn `R,C → 0` arrows are reversed in the new local convention. | `e_R=e_C=q/C`, `f_R+f_C=0`, `q̇=f_C` ⇒ `q̇=−q/(CR)`. |
| DC motor | Electrical 1-port signs `[+V,−R,−L,−GY.p1]`; mechanical `[+τ,−c,−J,+GY.p2]`; connect each pair; two-port laws change as above. The code's `+τ` is retained despite the drawn `jm → τ` arrow. | `L i̇=V−Ri−gω`, `Jω̇=τ+gi−cω` (`research/03-junctions-as-native-connections.md:95-100`). |
| 1-DoF | Replace today's `[+m,+d,+s,−f]` with `[−m,−d,−s,+f]`, the orientation drawn in TikZ. | `m v̇=F−cv−kq`, `q̇=v`; multiplying the old 1-junction sum by `−1` explains why the old signs worked. |
| 2-DoF | Use `m1_j1: [−m,−d,−s], out=−1`; `m2_j1: [−m,+f], out=+1`; `m2_sd: [−s,−d], out=+1`, then one 0-`connect` of their free ports (`research/03-junctions-as-native-connections.md:102-108`). Keep code placement, or deliberately revise the example and its goldens to match its drawing. | `q̇_b=v₁−v₂`, `m₁v̇₁=−k₁q₁−c₁v₁−B`, `m₂v̇₂=F+B`, as saved. |
| N-DoF | Keep each `sd_i` bonded **into** `0_i`: its 1-junction free port has `out=−1`; preceding mass junction has `out=−1` toward `0_i`; following mass junction receives from `0_i` with `out=+1`. All passive element ports are `−1`; connect the three free ports at each 0 node. Do not flip `sd_i → 0_i` to `0_i → sd_i` without changing state polarity. | `q̇_i=v_i−v_{i−1}` and the nearest-neighbor `−B_i+B_{i+1}` rows above. The prototype notes that reversing this branch gives mirrored velocities (`.scratch/mtk-upgrade/issues/04-prototype-connection-mechanism.md:23-25`). |
| Duffing | `[−m,−d,−s,−s3,+f]` already has the desired into-element orientation. | `m v̇=F−cv−kq−k₃q₃³`, `q̇=q̇₃=v`; redundant states remain. |

This algebra reproduces all six **code-built** equations. Separately, the accepted prototype on MTK 11.45.4 reports all six examples (nine golden cases) within maximum absolute simulation difference `2.3e−10` (`.scratch/mtk-upgrade/issues/04-prototype-connection-mechanism.md:11-25`, `docs/adr/0001-junctions-on-native-mtk-connections.md:10-22`). That is inherited evidence, not a fresh run in this analysis. The 2-DoF and DC TikZ differences remain documentation/model-choice questions; a solver match cannot resolve which drawing was intended.

## 6. Convention decision for v0.3

**Today, effectively:** `connect(a,b)` says `a → b`; the junction sum sees `+` on arrival and `−` on departure, while the shared variable is equated unsigned. A reversed junction argument rewrites `connect` order. This is consistent with the **junction** convention in the [20-sim tutorial](https://www.20sim.com/webhelp/modeling_tutorial_bond_graphs_frombondgraphtoequations.php) and with the DC motor output. It is **not a complete, uniformly enforced power-port convention**: one-port laws are written once regardless of arrow direction, the RC and 1-DoF drawings differ from local code polarity, and the two-port law assumes a particular pin-to-pout direction. Full energy accounting for arbitrary reversed element bonds is therefore uncertain without an explicit per-element flow convention and tests.

| Candidate | Consequence |
|---|---|
| **A. Accepted hybrid:** element bonds point into elements; flow positive into each component; native 0-`connect`; signed 1-junction ports; unsigned shared efforts/flows; GY/TF with one sign adjusted for two incoming ports. | Reproduces all code-built examples and golden trajectories in the reported MTK 11 prototype. `connect` order ceases to encode direction; 0-junction argument signs should be rejected or warned on; free 1-junction ports carry explicit orientation. Fits ADR 0001 and needs no MTK internals (`docs/adr/0001-junctions-on-native-mtk-connections.md:10-22`). |
| **B. Explicit components for both junctions:** same into-element ports and signed 1-junctions, but each 0-junction has its own component and ports. | Same algebra and no MTK internals; named 0 nodes are easier to inspect, at the cost of more alias equations and a guard against accidentally making an extra 0 node with a multi-member `connect` (`research/03-junctions-as-native-connections.md:114-125`). |
| **C. Fully native dual port for 1-junctions:** `PowerDual` swaps which variable is MTK `Flow`. | DC motor reproduced in research, but each element needs a port kind, sources change sign on dual ports, and every 0–1 bond needs an oriented adapter. Mixing port kinds can silently pair wrong variables (`research/03-junctions-as-native-connections.md:127-139`). The N-DoF chain is especially adapter-heavy; this is not the accepted design. |

**Decision evidence:** keep arrow signs in junction sums and keep junction equalities unsigned. Adopt A's explicit local port and 1-junction orientation rules for v0.3, as already accepted; correct or explicitly annotate the two contradictory drawings before using them as specification fixtures. Uncertain: exact internal route of the N-DoF nonlocal terms (the saved cancellations are certain), behavior of arbitrary reversed passive-element bonds in v0.2, and transformer sign behavior beyond the symbolic power-balance check—none of the six regression examples uses mTF (`TODO.md:27-45`).

STATUS: DONE
