# Model from examples/mechanical_msd_ndof.jl with the example's own N = 4 DoF.
function build_mechanical_msd_ndof(n = 4)
    @named m = Mass(m = 1.0)
    @named s = Spring(k = 1.0)
    @named d = Damper(c = 1.0)
    @named m1 = Junction1([-1, m], [-1, d], [-1, s])
    @named sd = Junction1([-1, s], [-1, d])
    M, cons = [], []
    for i = 2:n
        j0 = Junction0(sd; name = Symbol("m" * string(i) * "_j0"))
        j1 = Junction1([-1, m]; name = Symbol("m" * string(i) * "_j1"))
        push!(M, [j0, j1])
        push!(cons, connect(j0.power, j1.power))
    end
    push!(cons, [connect(M[i][2].power, M[i+1][1].power) for i = 1:(length(M)-1)]...)
    push!(cons, connect(m1.power, M[1][1].power))
    @named mdl = ODESystem(cons, t)
    mdl = compose(mdl, m1, vcat(M...)...)
    @named sys = simplifysys(mdl)
    prob = ODEProblem(sys, [M[1][2].m.power.f => 0.0, M[1][1].sd.s.q => 1.0], (0, 10))
    return sys, ["mechanical_msd_ndof" => prob]
end

check_example("mechanical_msd_ndof", build_mechanical_msd_ndof)
