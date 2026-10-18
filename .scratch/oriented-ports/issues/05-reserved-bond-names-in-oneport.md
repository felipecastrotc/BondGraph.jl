# 05: `e` and `f` are reserved inside `@oneport` bodies

**What to build:** since ADR 0002, `@oneport` defines the bond variables `e` and `f` before the body runs. A keyword or parameter named `f` (a frequency) or `e` clashes with them. MTK reports `ArgumentError: Found variables with duplicate names! f and f(t) have the same name but are not identical.`, which doesn't say why. At least document the reserved names in the `@oneport` docstring; possibly detect the clash in the macro and raise an error naming it.

**Blocked by:** None.

**Status:** needs-triage

## Reproduction

```julia
@oneport function SineSe(; name, A = 1.0, f = 1.0)
    @parameters A = A f = f
    [e ~ A * sin(2π * f * t)]
end
@named src = SineSe(A = 2.0, f = 0.5)   # ArgumentError: duplicate names f and f(t)
```

Script: `.scratch/api-usability/recheck/oneport_f_clash.jl`. Before ADR 0002 the bond variables were `power.e`/`power.f`, so the clash is new (not run on the old API).

## Acceptance criteria

- [ ] The `@oneport` docstring (and the user guide's element paragraph) names `e`, `f`, `power`, `pout` as reserved.
- [ ] Optional: the macro rejects a keyword or `@parameters`/`@variables` named `e` or `f` with a message naming it.

## Comments
