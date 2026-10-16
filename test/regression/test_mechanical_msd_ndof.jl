# Four masses; each coupling spring-damper bond points from the 0-junction to sd.
let n = 4
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0)
    @named d = Resistance(R = 1.0)
    @named m1 = Junction1([-1, m], [-1, d], [-1, s], [-1, :out])
    masses = [Junction1(:in, [-1, m], (i < n ? [[-1, :out]] : [])...;
                        name = Symbol(:m, i)) for i in 2:n]
    sds = [Junction1([-1, Compliance(C = 1.0, q = i == 2 ? -1.0 : 0.0, name = :s)],
                     [-1, d], :j0; name = Symbol(:sd, i)) for i in 2:n]
    prev = [m1; masses[1:end-1]]
    cons = [connect(prev[k].out, masses[k].in, sds[k].j0) for k in 1:n-1]
    @named model = System(cons, t; systems = [m1; masses; sds])
    cols = ["m1₊m₊power₊f(t)" => model.m1.m.power.f,
            "m1₊s₊q(t)" => model.m1.s.q]
    for i in 2:n
        push!(cols, "m$(i)_j1₊m₊power₊f(t)" => getproperty(model, Symbol(:m, i)).m.power.f)
        push!(cols, "m$(i)_j0₊sd₊s₊q(t)" => getproperty(model, Symbol(:sd, i)).s.q)
    end
    check_case("mechanical_msd_ndof", model, cols;
               flip = ["m$(i)_j0₊sd₊s₊q(t)" for i in 2:n])
end
