# CTMS: DC Motor Speed (System Modeling, PID Controller Design)

*Control Tutorials for MATLAB and Simulink* (CTMS), hosted by the University of Michigan, two pages of the Motor Speed example. Licensed CC BY-SA 4.0 (page footer).

- *DC Motor Speed: System Modeling*, <https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling>
- *DC Motor Speed: PID Controller Design*, <https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=ControlPID>

Read: the live site answered with a bot check (HTTP 403) on 2026-10-17, so both pages were read as Internet Archive snapshots: [System Modeling, 2022-11-19](https://web.archive.org/web/20221119140057/https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling) and [PID Controller Design, 2023-06-10](https://web.archive.org/web/20230610115150/https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=ControlPID) (both "Published with MATLAB 9.2"). Text, equations and MATLAB code were read; the figures (circuit sketch, block diagram, step-response plots) were not, so no value read off a plot is recorded here. Equation numbers are the pages' own.

Used by: `examples/control_dc_motor_pi.jl`. The same System Modeling page is linked from `examples/electrical_dc_motor.jl`.

## System Modeling page

- Setup: the input is the armature voltage V, the output the shaft speed θ̇. Rotor and shaft rigid; viscous friction, torque proportional to shaft speed.
- Parameters as stated: J = 0.01 kg·m² (rotor inertia), b = 0.1 N·m·s (viscous friction), Ke = 0.01 V/rad/s (EMF constant), Kt = 0.01 N·m/A (torque constant), R = 1 Ω, L = 0.5 H.
- Constant magnetic field ("armature-controlled motor"): T = Kt i (Eq. 1), back EMF e = Ke θ̇ (Eq. 2). In SI units Kt = Ke, and the page uses one K for both.
- Newton and Kirchhoff: J θ̈ + b θ̇ = K i (Eq. 3), L di/dt + R i = V − K θ̇ (Eq. 4). No load torque appears.
- Laplace: s(Js + b)Θ(s) = K I(s) (Eq. 5), (Ls + R)I(s) = V(s) − K s Θ(s) (Eq. 6). Eliminating I(s): P(s) = Θ̇(s)/V(s) = K / ((Js + b)(Ls + R) + K²), in (rad/s)/V (Eq. 7). MATLAB prints it as 0.01 / (0.005 s² + 0.06 s + 0.1001).
- State space with states (θ̇, i), input V, output θ̇ (Eqs. 8–9); MATLAB prints A = [−10 1; −0.02 −2], B = [0; 2], C = [1 0], D = 0.
- Design requirements: the uncompensated motor turns at 0.1 rad/s in steady state for 1 V (the page points to its System Analysis page for the simulation, not read). For a unit step in commanded speed: settling time < 2 s, overshoot < 5 %, steady-state error < 1 %.

## PID Controller Design page

- Restates Eqs. 5–7 as Eqs. 1–3 and the same parameters and requirements, "for a 1-rad/sec step reference". The loop is unity feedback: `sys_cl = feedback(C*P_motor,1)`. The block-diagram figure was not read.
- Controller: C(s) = Kp + Ki/s + Kd s = (Kd s² + Kp s + Ki)/s (Eq. 4).
- Sequence of designs, as stated:
  1. P, Kp = 100 (step plotted over 0–5 s): steady-state error and overshoot both too large; experimenting with Kp in Control System Designer confirms a proportional controller can't meet the requirements.
  2. PID, Kp = 75, Ki = 1, Kd = 1 (0–200 s): the steady-state error goes to zero, but reaching steady state takes far longer than 2 s, which the page attributes to the small Ki.
  3. PID, Kp = 100, Ki = 200, Kd = 1 (0–4 s): the error is removed much faster, but overshoot is greatly increased.
  4. PID, Kp = 100, Ki = 200, Kd = 10 (0–4 s): the larger Kd reduces the overshoot, and the page states that all the design requirements are met.
- The rise time, overshoot and settling time of each run are given only as plot annotations (not read). The page gives no PI design (Kd = 0).

## Relevance to BondGraphToolkit

- **Inference:** CTMS's motor parameters are exactly those of `examples/electrical_dc_motor.jl` (L = 0.5, R = 1, J = 0.01, c = b = 0.1, k = K = 0.01). That example adds what CTMS doesn't have: a 12 V constant source and a 1 N·m applied torque.
- **Inference:** Eqs. 3–4 are the bond graph of Karnopp Fig. 2.3d (see `karnopp-ch1-2.md`) with a gyrator of modulus K; Kt = Ke is the gyrator's power conservation (Karnopp Eq. 3.11).
- **Inference:** the denominator of Eq. 7 is 0.005 s² + 0.06 s + 0.1001, with poles at about −2.00 and −10.00. The control example uses a PI controller (C(s) with Kd = 0), with gains chosen by the example, not by CTMS: Kp = 20, Ki = 40 puts the PI zero at −2, near the slow pole, leaving a closed loop s² + 10 s + 2Kp (ζ ≈ 0.79). Simulated: 1.7 % overshoot, 2 % settling time 0.59 s, which meets the three requirements.
- **Inference:** CTMS's final PID uses an ideal derivative on the error, so a step reference gives an impulse in V; that is fine for `feedback`/`step` on transfer functions but would need a filtered derivative (e.g. MTKSL `Blocks.PID` with its `Nd`) in a time-domain model. This is why the example stays with PI.
