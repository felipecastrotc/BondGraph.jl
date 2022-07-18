using LinearAlgebra
using Statistics
using Symbolics
using SymPy
using Plots
using ProgressMeter
using LaTeXStrings
using SHA
using FLoops
using CodecZlib

pow(x, n) = prod([x for i in 1:n])

# ______________________________________________________________________________
# Viscosity models

# Viscosity models for water
# https://www.nist.gov/system/files/documents/srd/jpcrd382009101p.pdf
μₐ(T) = 2.414*1e-5*10^(247.8/(273.15 + T - 140));

# _______________________________________________________________________________
# Generate bases functions

# Based on SymPy - Faster solution with symbolic
# 1 → 33.821873 seconds (3.66 M allocations: 124.934 MiB, 0.19% gc time)
# 2 → 271.500200 seconds (22.37 M allocations: 600.757 MiB, 0.14% gc time) Eₚ ω q ρ D
# 3 → 3160.956391 seconds (227.74 M allocations: 6.073 GiB, 0.16% gc time) Eₚ Eₜ ω q ρ D
# function genbases(X, exps, var_units, out_unit, check_unit=true)
    
#     out_unit = string(out_unit);
#     # Generate symbolic variables
#     S = collect(symbols("x0:" * string(size(X, 1)), real=true));
#     units = Dict(zip(S, var_units));    # Translate variable to unit

#     B = [];     # Bases
#     Bₛ = [];    # symbolic bases
#     C = Iterators.product([exps for i in 1:size(X, 1)]...);
#     @showprogress 1 "Generating bases..." for i ∈ C
#     # p = Progress(length(C));
#     # Threads.@threads for i ∈ collect(C)
#         # Generate expression
#         expr = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(S)]);
#         #
#         if check_unit
#             # Get expression unit
#             expr_unit = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(var_units)]);
#             # expr_unit = subs(expr, units);
#             store = out_unit == string(expr_unit)
#         else
#             store = true
#         end

#         if store
#             push!(Bₛ, expr);    # Add the expression to the array
#             # Calculate the expression for the data
#             aux = prod(hcat([v.^i[j] for (j, v) in enumerate(eachrow(X))]...), dims=2);
#             push!(B, aux);
#         end
#         # next!(p)
#     end

#     return hcat(B...), Bₛ;
# end

# Based on Symbolocs library - way faster than sympy solution
# 854.341174 seconds (3.04 G allocations: 136.307 GiB, 5.41% gc time, 1.27% compilation time) -ΔP μ ω q ρ D [-6:6]
function genbases(var_dict, exps, out_unit, check_unit=true)
    # For adimensional the out_unit must be a string "a"

    # out_unit = string(out_unit);
    # Get settings from var_dict
    n_var = length(var_dict)
    x = keys(var_dict);
    units = values(var_dict);

    Bₛ = [];    # symbolic bases
    C = Iterators.product([exps for i in 1:n_var]...);
    @showprogress 1 "Generating bases..." for i ∈ C
    # p = Progress(length(C));
    # Threads.@threads for i ∈ collect(C)
        # Generate expression
        # expr = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(x)]);
        expr = prod([i[j] < 0 ? 1/pow(v,abs(i[j])) : pow(v,i[j]) for (j, v) in enumerate(x)])
        #
        if check_unit
            # Get expression unit
            # expr_unit = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(units)]);
            expr_unit = prod([i[j] < 0 ? 1/pow(v,abs(i[j])) : pow(v,i[j]) for (j, v) in enumerate(units)]);
            # expr_unit = subs(expr, units);
            if out_unit != "a"
                store = Symbolics.isequal(expr_unit, out_unit);
            else
                # Handling adimensional
                store = length(Symbolics.get_variables(expr_unit)) == 0;
            end
        else
            store = true;
        end

        if store
            push!(Bₛ, expr);    # Add the expression to the array
        end
        # next!(p)
    end

    return Bₛ;

end

function genbasesprl(var_dict, exps, out_unit, check_unit=true)
    # For adimensional the out_unit must be a string "a"

    # out_unit = string(out_unit);
    # Get settings from var_dict
    n_var = length(var_dict)
    x = keys(var_dict);
    units = values(var_dict);

    Bₛ = [];    # symbolic bases
    C = Iterators.product([exps for i in 1:n_var]...);
    # @showprogress 1 "Generating bases..." for i ∈ C
    # p = Progress(length(C));

    # Threads.@threads for i ∈ collect(C)
    @floop for i ∈ C
        # Generate expression
        # expr = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(x)]);
        expr = prod([i[j] < 0 ? 1/pow(v,abs(i[j])) : pow(v,i[j]) for (j, v) in enumerate(x)])
        #
        if check_unit
            # Get expression unit
            # expr_unit = prod([i[j] < 0 ? 1/v^(abs(i[j])) : v^i[j] for (j, v) in enumerate(units)]);
            expr_unit = prod([i[j] < 0 ? 1/pow(v,abs(i[j])) : pow(v,i[j]) for (j, v) in enumerate(units)]);
            # expr_unit = subs(expr, units);
            if out_unit != "a"
                store = Symbolics.isequal(expr_unit, out_unit);
            else
                # Handling adimensional
                store = length(Symbolics.get_variables(expr_unit)) == 0;
            end
        else
            store = true;
        end

        if store
            @reduce(Bₛ = vcat(Num[], [expr]));
        end
        # next!(p)
    end

    return Bₛ;

end

# Generate bases functions wihtout symbolics
function genbases(X, exps)
    # Gen bases for dimensionless problems
    B = [];     # Bases
    C = Iterators.product([exps for i in 1:size(X, 1)]...);
    for i ∈ C
        aux = prod(hcat([v.^i[j] for (j, v) in enumerate(eachrow(X))]...), dims=2);
        push!(B, aux);
    end

    return hcat(B...);
end

# Function to regenerate bases according to the dataframe passed and symbolic
# array of symbolic expressions 
function regenbases(Bₘₛ, dff, s2i, var_array::Matrix{SymPy.Sym})
    # Empty bases
    Bₙ = [];
    @showprogress 1 "Generating bases..." for expr ∈ Bₘₛ
        if length(collect(expr.free_symbols)) != 0
            idx = sort([s2i[s] for s in collect(expr.free_symbols)]);
            # Get dataframe column names
            cols = [string(var_array[i]) for i in idx];
            # Get data from the dataframe
            dt = [convert(Array, dff[!, c]) for c in cols];
            # Check if it is not a 
            # Get function from SymPy
            foo = lambdify(expr);
            # Store data
            out = foo.(dt...);
        else
            # In case of no arguments are used e.g. x^0*y^0*z^0
            out = ones(size(dff, 1));
        end
        push!(Bₙ, out);
    end
    return hcat(Bₙ...);
end

function regenbases(Bₘₛ, dff, var_array::Matrix{SymPy.Sym})
    # Model symbols
    S = collect(symbols("x0:" * string(length(var_array)), real=true));

    # Translate symbol to index
    s2i = Dict(zip(S, 1:length(S)));

    return regenbases(Bₘₛ, dff, s2i, var_array)

end 

function regenbases(Bₛ, dff, var_array::Matrix{Num})

    # Get data
    b = dff[!, vec([string(i) for i in var_array])] |> Array;

    B = [];
    # Create variables symbols
    vₛ = [Symbol(v) for v in var_array];

    # Empty bases
    
    p = Progress(length(Bₛ));
    for bₛ in Bₛ
        # Convert symbolic to expression
        expr = Symbolics.toexpr(bₛ);
        # Create function from expression
        foo = quote
            $([:($v = x[$i]) for (i, v) in enumerate(vₛ)]...) 
            return $(expr)
        end
        f = @eval (x) -> $foo ;
        # Evaluate function
        # push!(B, f.(eachrow(b)));
        push!(B, Base.invokelatest.(f, eachrow(b)));
        next!(p);
    end

    return hcat(B...);

end 

function regenbases(Bₛ, dff, var_array, cache::String)
    
    # Get data
    b = dff[!, vec([string(i) for i in var_array])] |> Array;
    # Check cache
    if length(cache) > 1
        BB = getbases(cache, b);
    end
    if length(BB) == 0
        
        BB = regenbases(Bₛ, dff, var_array::Matrix{Num})

        # Store cache
        if length(cache) > 1
            updatecache(cache, b, BB);
        end
    end

    return BB;
end

function regenbases(Bₘₛ, dff)
    # Empty bases
    Bₙ = [];
    @showprogress 1 "Generating bases..." for expr ∈ Bₘₛ
        # Get dataframe column names
        # cols = sort([string(s) for s in collect(expr.free_symbols)]);
        cols = [string(s) for s in collect(expr.free_symbols)];
        # Get data from the dataframe
        dt = [convert(Array, dff[c]') for c in cols];
        # Create ordered symbols from the free symbols
        S = collect(symbols("x0:" * string(size(collect(expr.free_symbols), 1)), real=true));
        # Translate expression symbols to ordered symbols
        o = Dict( e => s for (e, s) in zip(collect(expr.free_symbols),S));
        expr = subs(expr, o);
        # Get function from SymPy
        foo = lambdify(expr);
        # Store data
        push!(Bₙ, vec(foo.(dt...)));
    end
    return hcat(Bₙ...);
end

# _______________________________________________________________________________
# Trim model

# Simplify the model basing on the relative contribution to the final results,
# it fits a OLS at each coefficient removal.
function trim_model(B, Bₛ, ỹ, trim=0.02)

    # Least square fit
    θ = B\ỹ;
    
    # Create a temporary base array to be trimmed
    Bₜ = B;
    Bᵢ = collect(1:size(B, 2));     # Index tracking
    for i ∈ 1:10
        # Analyse contribution of each parameter
        Ψ = (Bₜ'.*θ)';      # Matrix with the indivudual contribution
        ψ = abs.(Ψ)./sum(abs.(Ψ), dims=2);     # Relative contribution
        # Remove Inf and NaN cases
        replace!(ψ, Inf=>0);
        replace!(ψ, NaN=>0);
        
        # Weighted mean based on error
        err = abs.(ỹ .- sum(Ψ, dims=2));
        err /= sum(err);

        # # Get all variables contribution
        ψₒ = vec(sum(ψ.*err, dims=1));
        ψᵢ = sortperm(ψₒ);      # Indexes sorted

        # # Get mask
        cdf = cumsum(ψₒ[ψᵢ]);
        κ = ψᵢ[cdf .> trim];

        if sum(cdf .< trim) == 0
            break
        end

        # Apply mask to the bases
        Bᵢ = Bᵢ[κ];
        Bₜ = Bₜ[:, κ];
        θ = Bₜ\ỹ;
    end
    
    θ = Bₜ\ỹ;

    return return Bₜ, Bₛ[Bᵢ], θ;
end

# Simplify the model coefficients basing on the CDF
function trim_model(θ, trim=0.001)
    
    # Sort the model coefficients
    i = sortperm(abs.(θ));

    # Calculate the cummulatinve distribution function
    cdf = cumsum(abs.(θ[i]))/sum(abs.(θ));

    # Get indexes that are higher than the threshold
    i = i[findmax(diff(cdf .> trim))[2]+1:end];

    return i
end

# _______________________________________________________________________________
# Model symbolics

function get_equations(θ, Bₛ)
    
    # Model symbols
    bₛ = sympify.([string(b) for b in Bₛ]);
    
    # Model constants
    Θₖ = collect(symbols("k1:" * string(length(bₛ) + 1), real=true));
    Θₖ = reshape(Θₖ, (1, length(Θₖ)));

    # Generate equation with model constants
    eqk = (Θₖ*bₛ)[1];

    # Generate equation with model coefficients
    θ = reshape(θ, (1, length(θ)));
    eqc = (θ*bₛ)[1];

    return eqc, eqk;

end

function group_expr(expr, K)

    # Set of constants
    # K = [k₁, k₂, k₃, k₄, k₅, k₆, k₇, k₈, k₉, k₁₀, k₁₁];
    # Remove constant values
    S = Dict(k => 1 for k in K);
    expr = SymPy.expand(expr);
    atoms = [];
    # Get atoms wihtout constants and multiplying numbers
    for a in expr.args
        if length(a.args) > length(a.free_symbols) 
            a = subs(a, a.args[1] => 1);
        end
        push!(atoms, subs(a, S));
    end
    # Remove redundant atoms
    atoms = collect(Set(atoms));
    # Create a new symbol for each atom
    E = collect(symbols("x0:" * string(length(atoms)), real=true));
    # Assign each atom with a symbol
    C = Dict(a => e for (a, e) in zip(atoms, E));
    expr = subs(SymPy.expand(expr), C);
    # Collect atoms
    for (a, e) in zip(atoms, E)
        expr = collect(expr, e => a);
    end
    # Substitute symbols
    C = Dict(e => a for (a, e) in zip(atoms, E));
    return subs(expr, C)
end

# _______________________________________________________________________________
# Load and store utils

# Generic functions
function store(name, data)
    filename = "./models/" * name * ".out";
    open(filename, "w") do f
        serialize(f, data);
    end;
end

function load(name)
    filename = "./models/" * name * ".out";
    if isfile(filename)
        data = open(filename) do f
            return deserialize(f);
        end;
    else
        data = Dict();
    end
    return data;
end

# Cache functions
function getcached(name)
    return load(name*"_cached")
end

function storecache(name, data)
    return store(name*"_cached", data)
end

# Model functions
function storemodel(name, Bₛ, var_dict)
    data = Dict("bases" => Bₛ, "vars" => var_dict)
    store(name, data);
end

function loadmodel(name)
    data = load(name);
    return data["bases"], data["vars"]
end

# _______________________________________________________________________________
# Cache
function df2hash(dt)
    return bytes2hex.(sha256.(eachrow(convert(Array, dt)) .|> string));
end

function getbases(name, dt)
    
    cached = getcached(name);
    cached_hash = keys(cached);

    hashes = df2hash(dt);

    if sum([h in cached_hash for h in hashes]) == length(hashes)
        B = hcat([cached[h] for h in hashes]...)';
    else
        B = [];
    end

    return B
end

function updatecache(name, dt, B)

    cached = getcached(name);
    hashes = df2hash(dt);

    for (h, b) in zip(hashes, eachrow(B))
        cached[h] = b;
    end

    storecache(name, cached);
end

# _______________________________________________________________________________
# Complexity
function compressratio(sym)
    if typeof(ω/μ) == Num
        expr = Symbolics.toexpr(sym) |> string
    else
        expr = sym |> string
    end

    exp_bin = Vector{UInt8}(expr);
    exp_zlib = transcode(ZlibCompressor, expr);

    bexpr = join([bitstring(hton(i)) for i in exp_bin])
    bzlib = join([bitstring(hton(i)) for i in exp_zlib])

    bolexpr = [parse(Bool, i) for i in bexpr];
    bolzlib = [parse(Bool, i) for i in bzlib];

    return length(bolzlib)/length(bolexpr);
end

function compressratioLZW(expr)
    # if typeof(ω/μ) == Num
    #     expr = Symbolics.toexpr(sym) |> string
    # else
    #     expr = sym |> string
    # end

    # length(expr)
    # exp_bin = Vector{UInt8}(expr);
    exp_lzw = compressLZW(expr);

    # bexpr = join([bitstring(hton(i)) for i in exp_bin])
    # blzw = join([bitstring(hton(i)) for i in exp_lzw])

    # bolexpr = [parse(Bool, i) for i in bexpr];
    # bollzw = [parse(Bool, i) for i in blzw];

    # return length(bollzw)/length(bolexpr);
    # return length(exp_lzw)/length(exp_bin);
    return length(exp_lzw)/length(expr);
end

function expr2ascii(Bₛ)
    
    # Alphabet to be used
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVXYZW";
    alphabet = alphabet*lowercase(alphabet);

    # Get unique variables used in all expressions
    vars = unique(reduce(vcat, Symbolics.get_variables.(Bₛ)));
    
    # Remove the variables already being used  from the alphabet
    str_un = [occursin(i |> string, alphabet) for i in vars];
    str_av = setdiff(Set(alphabet), Set([ i |> string for i in vars[str_un]]));
    str_arr = collect(str_av);      # Available chars in an array format

    # Get variables needing to be replaced
    msk = .!str_un;
    # Select variables to replace
    vars_replace = Symbolics.variable.(str_arr[1:sum(msk)]);
    # Create a dictionary to replace the variables
    vars_dct = Dict(i => j for (i, j) in zip(vars[msk], vars_replace));

    # Replace expressions
    return [substitute(Bₛ[i], vars_dct) for i in 1:size(Bₛ, 1)];
end

function compressLZW(decompressed::String)
    dictsize = 256
    dict     = Dict{String,Int}(string(Char(i)) => i for i in 0:dictsize)
    result   = Vector{Int}(undef, 0)
    w        = ""
    for c in decompressed
        wc = string(w, c)
        if haskey(dict, wc)
            w = wc
        else
            push!(result, dict[w])
            dict[wc]  = dictsize
            dictsize += 1
            w        = string(c)
        end
    end
    if !isempty(w) result[1] = dict[w] end
    return result
end

function compressratioLZWd(expr::String)

    dictsize = 256
    dict = Dict{String,Int}(string(Char(i)) => i for i in 0:dictsize)
    w = ""
    i = 1
    for c in expr
        wc = string(w, c)
        if haskey(dict, wc)
            w = wc
        else
            i += 1;
            dict[wc]  = dictsize
            dictsize += 1
            w = string(c)
        end
    end
    if !isempty(w) i -= 1 end
    return i/length(expr)
end

function compressratioLZWd(d::Vector{Bool})

    dict = Dict{Vector{Bool},Int}([true] => 1, [false] => 2)
    dictsize = 2;

    w = Vector{Bool}(undef, 0)
    i = 0;
    for c in d
        wc = vcat(w, c)
        if haskey(dict, wc)
            w = wc
        else
            i += 1
            dictsize += 1
            dict[wc] = dictsize
            w = [c];
        end
    end

    return (i)/length(d)
end


# _______________________________________________________________________________
# Utils

function plot_result(y, ŷ, df, title="", size=(800, 800))
    
    l = @layout [a ; b]

    # Plot ΔP approx vs ΔP truth 
    p1 = scatter(y, ŷ, xlabel="Model [Pa]", ylabel="Truth [Pa]", title=title, legend=false);

    # Plot ΔP x q
    p2 = scatter(df[!, "q"], y, xlabel=L"q \, [m^3/s]", ylabel=L"\Delta P \, [Pa]", label="Truth");
    scatter!(p2, df[!, "q"], ŷ, label="Model");

    ax = plot(p1, p2, size=size, layout = l);

    return ax;

end