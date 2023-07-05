
"""
    equalityeqs(con)

Create a list of equality equations from a given array of expressions.

# Arguments:
- con: Array of expressions.

# Returns:
- Array of equality equations.

# Example:
```julia-repl
julia> con = [a, b, c]

julia> eqs = equalityeqs(con)
julia> # eqs = [a ~ b, a ~ c]

julia> con = [a]

julia> eqs = equalityeqs(con)
julia> # eqs = []
```
"""
function equalityeqs(con)
    base = 0.0
    if length(con) > 1
        base += con[1]
    else
        return []  # Return an empty array if there's only one expression
    end

    eqs = Equation[]
    for i = 2:length(con)
        push!(eqs, base ~ con[i])  # Create an equality equation between base and each expression
    end

    return eqs
end

"""
    sumvar(con)

Create a sum variable equation from a given array of expressions.

# Arguments:
- con: Array of expressions.

# Returns:
- Sum variable equation if the array is non-empty, otherwise an empty array.

# Example:
```julia-repl
julia>  con = [a, b, c]

julia>  eq = sumvar(con)
# eq = 0 ~ (a + b + c)

julia> con = []

julia> eq = sumvar(con)
# eq = []
```
"""
function sumvar(con)
    if length(con) > 0
        return 0 ~ sum(con)  # Create a sum variable equation with 0 on the left-hand side and the sum of expressions on the right-hand side
    else
        return []  # Return an empty array if the input array is empty
    end
end

"""
    flatinput(ps)

Flatten and process input systems and signs.

# Arguments:
- ps: Array of input systems or input vectors.

# Returns:
- Tuple containing the flattened input systems and corresponding signs.

# Example:
```julia-repl
julia> ps = [sys1, [-1, sys2], sys3]

julia> subsys, signs = flatinput(ps)
# subsys = [sys1, sys2, sys3]
# signs = [1, -1, 1]
```
"""
function flatinput(ps)
    subsys = ModelingToolkit.AbstractSystem[]  # Array to store the flattened input systems
    signs = Int[]  # Array to store the corresponding signs

    for (i, p) in enumerate(ps)
        if isa(p, ModelingToolkit.AbstractSystem)
            push!(subsys, p)
            push!(signs, 1)  # Add a positive sign (1) to the signs array (default)
        elseif isa(p, AbstractVector)
            # Check if the AbstractVector has two elements where the first is an integer and the second is an AbstractSystem
            if (length(p) == 2) &&
               isa(p[1], Int) &&
               isa(p[2], ModelingToolkit.AbstractSystem)
                push!(subsys, p[2])
                push!(signs, p[1])  # Add the sign (integer) to the signs array
            else
                throw(
                    DomainError(
                        "The input number " *
                        string(i) *
                        " is not valid, a valid element has a two elements AbstractVector, where the first is an integer and the second is an AbstractSystem",
                    ),
                )
            end
        else
            throw(
                DomainError(
                    p,
                    "The input number " *
                    string(i) *
                    " is not valid. The valid options are: AbstractSystem or a two-element AbstractVector.",
                ),
            )
        end
    end

    return subsys, signs  # Return the flattened input systems and corresponding signs as a tuple
end

"""
    isindependent(var::Num)

Check if a variable is independent.

# Arguments:
- var::Num: The variable to check.

# Returns:
- true if the variable is independent.
- false if the variable is not independent.
"""
function isindependent(var::Num)
    if isvariable(var)
        !istree(unwrap(var))
    else
        true
    end
end
