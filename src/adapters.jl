# Adapters to MTK standard-library (MTKSL) networks. An adapter is a one-port element
# whose bond conventionally points into it: the MTKSL network on its MTKSL side is the
# element at the end of that bond, with `e` the effort across the network and `f` the
# flow into it. Like every element it has the ports `power` and `pout` (ADR 0002).
# MTKSL counts its through-variables positive into each component, so the adapter's own
# through-variable is negated. When the MTKSL side supplies power, the power along the
# bond is simply negative.
#
# The mapping is research 11's (`.scratch/mtk-upgrade/research/11-coupling-with-mtk-stdlib.md`):
# effort = across-variable in the impedance-analogy domains (electrical, hydraulic);
# flow = rate of the across-variable in the mobility-analogy ones (mechanical).
#
# The adapters are plain `@component`s, not `@oneport`s: `@oneport` only adds the
# `input`/`output` signal connectors, and an adapter also needs an MTKSL connector (pin,
# flange, port). They take their ports and bond variables from `bond_ports`. They are hand-written elements, as in "Writing an element without the
# macros" (`examples/elements_oneport.jl`).

"""
    ElectricalAdapter(; name)

Bond to an MTKSL electrical network connected between pins `p` and `n`:
`e = p.v - n.v` (voltage across the network) and `f` = current into the network at `p`.
Its bond conventionally points into the adapter; a grounded network uses
`connect(adapter.n, ground.g)`.
"""
@component function ElectricalAdapter(; name)
    power, pout, e, f, eqs = bond_ports()
    @named p = Pin()
    @named n = Pin()
    System([eqs; e ~ p.v - n.v; f ~ -p.i; 0 ~ p.i + n.i], t; name, systems = [power, pout, p, n])
end

"""
    TranslationalAdapter(; name)

Bond to an MTKSL `Translational` (velocity-based) network at `flange`: `e` = force on the
network, `f = flange.v`. Its bond conventionally points into the adapter.
"""
@component function TranslationalAdapter(; name)
    power, pout, e, f, eqs = bond_ports()
    @named flange = Translational.MechanicalPort()
    System([eqs; e ~ -flange.f; f ~ flange.v], t; name, systems = [power, pout, flange])
end

"""
    TranslationalPositionAdapter(; name, s = nothing)

Bond to an MTKSL `TranslationalPosition` network at `flange`: `e` = force on the network,
`f = ds/dt`. `s` is the flange's initial position; leave it `nothing` when the MTKSL side
already sets it. Its bond conventionally points into the adapter.
"""
@component function TranslationalPositionAdapter(; name, s = nothing)
    power, pout, e, f, eqs = bond_ports()
    @named flange = TranslationalPosition.Flange()
    @variables s(t) = s [guess = 0.0]
    System([eqs; s ~ flange.s; D(s) ~ f; e ~ -flange.f], t; name, systems = [power, pout, flange])
end

"""
    RotationalAdapter(; name, phi = nothing)

Bond to an MTKSL `Rotational` network at `flange` (fixed support): `e` = torque on the
network, `f = dφ/dt`. `phi` is the flange's initial angle; leave it `nothing` when the
MTKSL side already sets it. Its bond conventionally points into the adapter.
"""
@component function RotationalAdapter(; name, phi = nothing)
    power, pout, e, f, eqs = bond_ports()
    @named flange = Rotational.Flange()
    @variables phi(t) = phi [guess = 0.0]
    System([eqs; phi ~ flange.phi; D(phi) ~ f; e ~ -flange.tau], t; name, systems = [power, pout, flange])
end

"""
    HydraulicAdapter(; name)

Bond to an MTKSL `IsothermalCompressible` network at `port`: `e = port.p` and `f` = volume
flow into the network at the fluid's reference density, `-port.dm / port.ρ`. Connect the
`HydraulicFluid` in its own statement (`connect(fluid, src.port)`), not in the same
`connect` as the adapter. Its bond conventionally points into the adapter.
"""
@component function HydraulicAdapter(; name)
    power, pout, e, f, eqs = bond_ports()
    @named port = IsothermalCompressible.HydraulicPort()
    System([eqs; e ~ port.p; f ~ -port.dm / port.ρ], t; name, systems = [power, pout, port])
end
