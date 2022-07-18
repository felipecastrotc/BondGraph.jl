using CSV
using LinearAlgebra
using Plots
using LaTeXStrings
using ScikitLearn
using Statistics
using Serialization
using Random
using Symbolics
using LsqFit

include("lib_sparse_units.jl")


# ______________________________________________________________________________
# Prepare data
# Get symbolic variables
@variables ρ, μ, L, V, ΔP, ϵ, a, A, d, Q
# Units
@variables m, kg, s, K
units = Dict(ρ => kg/m^3, μ => kg/(m*s), L => m, V => m/s, ΔP => kg/(m*s^2), ϵ => m, a => m/s, A => m^2, d => m, Q => m^3/s)

# var_array = [ρ, μ, d, L, ϵ, a, A, V];
var_array = [ρ, μ, d, L, ϵ, A, V];

# Create an array of units
var_units = [units[s] for s in var_array];
var_dict = Dict(zip(var_array, var_units));

exps = collect(-5:5);

# out_unit = (units[μ]*units[L]/(units[d]^4))
out_unit = (units[ρ]*units[L]/units[d])*units[V]
# @time Bₛ = genbases(var_dict, exps, out_unit, true);      # HEAVY!!!!
@time Bₛ = genbasesprl(var_dict, exps, out_unit, true);      # HEAVY!!!!

Bₛ
filename = "./data/model_damper.out"
# open(filename, "w") do f
#     serialize(f, Bₛ);
# end;

Bₛ = open(filename, "r") do f
    return deserialize(f);
end;

using ProgressMeter


B = [];
# Create variables symbols
vₛ = [Symbol(v) for v in var_array];

B = []
p = Progress(length(Bₛ));
# Empty bases
for bₛ in Bₛ
    # Convert symbolic to expression
    expr = Symbolics.toexpr(bₛ);
    # Create code from expression
    tmp = quote
        b = zeros(size(x, 2))
        for (i, xc) in enumerate(eachcol(x))
            $([:(local $v) for (i, v) in enumerate(vₛ)]...)
            $([:($v = xc[$i]) for (i, v) in enumerate(vₛ)]...)
            b[i] = $(expr)
        end
        push!(B, b)
    end
    eval(tmp)
    next!(p);
end

B = hcat(B...)

idx = 1:length(Bₛ)
θ = B\dlist
θ = B\yd

Bₛ[3237]

idx[θ .> 1e-2]
θ[θ .> 1e-3]
