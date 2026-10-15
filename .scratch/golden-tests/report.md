# Golden-data regression tests (pre-MTK-upgrade safety net)

## What was built
- `test/regression/test_<example>.jl`: `build_<example>()` rebuilds the model exactly as `examples/<example>.jl` and returns `(sys, ["case" => ODEProblem, ...])`, then calls `check_example`.
- `test/regression/helpers.jl`: the only place that knows state naming (`state_names = string.(states(prob.f.sys))`), simulation (Tsit5, abstol=reltol=1e-10, retcode checked, `sol(TS)` with `TS = 0:0.5:10`, 21 points) and CSV read/write/compare (DelimitedFiles, columns matched by name, elementwise `isapprox` rtol=1e-6 atol=1e-9; a missing reference state fails naming it, extra states only warn).
- `test/regression/reference/<case>.csv` (header `t` + `string(state)`) and `<example>_equations.txt` (simplified equations before the forcing substitution; not asserted).
- `test/regression/generate_references.jl`; `Regression` testset in `test/runtests.jl`; DelimitedFiles in `[extras]`, `[targets].test`, `[compat] = "1"`.
- Stack: Julia 1.10.12, MTK 8.71.2, Symbolics 5.8.0, SymbolicUtils 1.4.0, OrdinaryDiffEq 6.58.0, DifferentialEquations 7.11.0.

## Cases (tspan (0,10), 21 samples)

| case | states |
|---|---|
| electrical_rc | 1 |
| electrical_dc_motor (bond-graph model only) | 2 |
| mechanical_duffing_unforced / _forced (F = 8cos(0.5t)) | 3 |
| mechanical_msd_1dof_unforced / _forced (F = sin t) | 2 |
| mechanical_msd_2dof_unforced / _forced (F = sin t) | 4 |
| mechanical_msd_ndof (N = 4, the example's own n) | 8 |

The examples use `solve(prob)`; the tests pin Tsit5. No example was skipped.

## Physics review (reported, not fixed)
- Energy decays in all damped cases; RC follows exp(-t); the 1-DoF forced case approaches -cos(t); DC motor steady state i≈11.888, ω≈11.189 with signs matching L di/dt = V − R i − gω and J dω/dt = τ + g i − cω.
- **N-DoF:** equations contain cancelling terms (m3/m4 spring pairs in the m1 equation, the m4 pair in m2). The 0-junction connection sets merge non-adjacent elements; physics is correct after cancellation, but this is the most fragile case under the upgrade.
- **Duffing:** redundant state (`s₊q` and `s3₊q` both integrate `m₊power₊f`), equal only because the initial conditions are. MTK may eliminate one after the upgrade.
- The `TODO.md` sign-convention question is unresolved.

## Commands

```sh
julia --project -e 'using Pkg; Pkg.test()'
julia --project -e 'using Pkg; Pkg.test(test_args = ["regenerate"])'   # only on a known-good stack
```

If MTK renames states, adapt `state_names` in `helpers.jl` (or map old to new names there) instead of regenerating.

## Test summary
`Pkg.test()` passes: base 20, utilities 41, ports and system 155, macros 5, Regression 38 (10 cases). Total 259. Verified independently by the parent session. A renamed CSV column fails with the expected missing-state message.

STATUS: DONE
