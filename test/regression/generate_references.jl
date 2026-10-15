# Regenerate test/regression/reference/ from the current stack:
#   julia --project -e 'using Pkg; Pkg.test(test_args = ["regenerate"])'
using BondGraph
include("helpers.jl")
REGENERATE[] = true
for f in filter(startswith("test_"), readdir(@__DIR__))
    include(f)
end
