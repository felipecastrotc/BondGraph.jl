# Model from examples/electrical_rc.jl: R and C on one 0-junction, bonds drawn into them.
let
    @named C = Compliance(C = 1.0, q = 1.0)
    @named R = Resistance(R = 1.0)
    @named model = Junction0([-1, R], [-1, C])
    check_case("electrical_rc", model, ["C₊q(t)" => model.C.q])
end
