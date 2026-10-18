# Spec: oriented ports

Status: ready-for-agent

Sources: [ADR 0002](../../docs/adr/0002-oriented-ports.md), ticket 12 (`.scratch/v0.3-impl/issues/12-manual-checks.md`), `GLOSSARY.md` (incoming/outgoing port).

## Problem

A plain `connect` can't make a bond out of an element: elements only have a port whose bond points into them, so `connect(se.power, tf.p1)` silently negates the source's flow and power (and an Sf's physics). Bond direction can only be given as a junction member's `[-1, x]`, which can't be used with `connect` or submodels.

## Solution

Every element has an incoming port `power` and an outgoing port `pout`, and the variables `e`, `f` of its bond along the half-arrow. A bond is a `connect` from an outgoing port (`pout`, `p2`, a junction's outgoing free port) to an incoming one (`power`, `p1`, a junction's incoming free port). Junction members: `x` (bond into x), `x => :pout` (bond out of x), `:a` (free port, bond leaving the junction), `:a => :in` (bond entering it). `[-1, x]` raises an error. `bondgraph` raises an error when an element is bonded through both ports.

## Decisions

- Port names `power` / `pout` (user, 2026-10-18, option (a)).
- Junction element member with its outgoing port: `x => :pout`. Found while implementing: a junction must own its elements as subsystems (`m.se`), and `x.pout` is a namespaced connector that doesn't lead back to `x`, so `Junction1(x.pout, …)` can't work. `x => :power` is accepted as the explicit form of `x`; `:a => :out` of `:a`.
- Error for both ports in `bondgraph` only, no compile wrapper (user).
- Two-ports keep `p1` / `p2` (user).
- Results are read from `x.e`, `x.f`; a port's raw `f` keeps MTK's sign.
- `@oneport` bodies use `e`, `f` (the bond's variables); `power.e` / `power.f` stay valid MTK variables but `power.f` is the bond flow only when the bond uses `power`.

## Out of scope

Signals in drawings (ticket 20), `bg_equations` observed (22), derivative-causality docs (21).

## Acceptance (whole feature)

- `Pkg.test()` green; the v0.2 regression data matches, not regenerated.
- `Se ⇀ TF ⇀ I` and `Sf ⇀ TF ⇀ C` with plain `connect`s give the closed-form answers (tests).
- Every example runs with `julia --project=docs`; the docs build.
