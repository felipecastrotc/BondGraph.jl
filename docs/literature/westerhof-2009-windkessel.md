# Westerhof, Lankhaar & Westerhof 2009: The arterial Windkessel

N. Westerhof, J.-W. Lankhaar and B. E. Westerhof, *The arterial Windkessel*, Med Biol Eng Comput 47: 131–141, 2009 (published online 10 June 2008; special issue, review). Open access (CC BY-NC), doi:10.1007/s11517-008-0359-2.

Read: the Springer HTML full text (link.springer.com/article/10.1007/s11517-008-0359-2), all of §1–§8, Figs. 3 and 5 viewed as images, Tables 1–2; the other figure captions read as text only. The paper is not in PMC (PMC2644372 is a different article). Section, figure and table numbers below are the journal's; the HTML has no page numbers, so none are given.

Used by: `examples/biomedical_windkessel.jl` (three-element Windkessel).

## Two-element (Frank) Windkessel (§1, §2)

- Frank formulated the two-element model: a resistance and a compliance. Weber had earlier compared the elastic large arteries with the air chamber (Windkessel) of a fire engine, which turns pulsatile pumping into a nearly steady outflow. (§1, Fig. 1)
- Peripheral resistance R: mostly in the smallest arteries and arterioles. R = (P_ao,mean − P_ven,mean)/CO ≈ P_ao,mean/CO, with CO the cardiac output. (§1, unnumbered equation)
- Total arterial compliance C = ΔV/ΔP, mainly from the large (conduit) arteries; hard to measure directly because volume leaks through the periphery. (§1)
- In diastole (valve closed) the model predicts an exponential pressure decay, P_dia(t) = P_es e^(−t/RC), P_es the end-systolic pressure; RC characterises the arteries, the heart period T the heart. (§2, Fig. 3)
- Its input impedance falls to negligible values with phase −90° at high frequencies, unlike the measured aortic impedance, and it predicts systolic pressure poorly. (§2, Fig. 4; Table 1: "total waveshape: poor")

## Three-element Windkessel (§2)

- Adds the characteristic impedance of the proximal aorta, Z_c, in series at the entrance; the measured impedance modulus levels off at Z_c at high frequencies. (§2, Fig. 4)
- Topology (Fig. 5, hydraulic and electrical forms, viewed): heart → Z_c → node with C to ground in parallel with R to the veins. The electrical drawing is Z_c in series, then C ∥ R.
- Z_c = PWV·ρ/A (pulse wave velocity, blood density, proximal aortic area); it links the lumped model to wave travel. (§2, Fig. 5 caption)
- Z_c has the dimension of a resistance and is usually drawn as one, but it is not a resistance: it only has meaning for oscillations. With a resistor for Z_c, mean pressure over mean flow is R + Z_c, not R, an error at low frequencies. (§2; Table 1)
- In the systemic circulation of all mammals, Z_c is about 5–7% of R, so that error is small. (§2, citing their ref. 72)
- Fed with measured aortic flow, it gives pressures close to measured ones, missing high-frequency detail (inflection point, augmentation). (§2, Fig. 6)
- A nonlinear (pressure-dependent C) version is not preferred over the linear one (§2, quoting Fogliardi et al.); "the three-element Windkessel suffices in most studies".
- Diastolic decay time is RC for both the 2- and 3-element models, when the fit starts about 10% of the heart period after valve closure. (§5, decay time method)
- Volume balance used by the two-area method: ∫Q dt = C (P₁ − P₂) + (1/R) ∫P dt over an interval. (§5)

## Four-element Windkessel and others (§3)

- Adds a total arterial inertance L in parallel with Z_c (Fig. 5) to fix the low-frequency error; L is hard to estimate. Others put an inertance in series with Z_c. (§3, Table 1)

## What the paper does not give

- **No numerical parameter values** (no R, C or Z_c in mmHg·s/mL or other units), only the ratio Z_c ≈ 5–7% of R. Tables 1–2 are qualitative or sensitivities; Fig. 3 has no axis values.
- **No written ODE** for the three-element model; only the diastolic decay, the input-impedance form (a + jωb)/(1 + jωc), stated in §3 as shared by the three-element model and a related model of Burattini, and the volume balance of §5.
- **No inflow waveform**: the examples use measured aortic flow (Fig. 6) or an isolated heart (Figs. 7–8).

## Relevance to BondGraphToolkit

- **Inference:** with pressure as effort and volume flow as flow, Fig. 5 is a bond graph Sf(Q) → 1 → {R:Z_c, 0 → {C, R:R}}: Z_c carries the whole inflow (1-junction), C and R share the pressure after Z_c (0-junction), and the veins are the pressure reference.
- **Inference:** the hand-derived equations from Fig. 5's circuit are C dP_c/dt = Q − P_c/R and P = P_c + Z_c Q, with P the entrance pressure and P_c the pressure across C. Q = 0 gives the paper's diastolic decay with time constant RC, and steady state gives mean P = (R + Z_c)·mean Q, as Table 1 states.
- Example choices, not from the paper: the half-sine systolic inflow, and the values R = 1.0 mmHg·s/mL, Z_c = 0.06 R (inside the paper's 5–7%), C = 1.5 mL/mmHg, T = 0.8 s, systole 0.3 s, SV = 70 mL.
