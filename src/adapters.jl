# Adapters to MTK standard-library (MTKSL) networks. An adapter is a one-port element
# whose bond conventionally points into it: the MTKSL network on its MTKSL side is the
# element at the end of that bond, with `power.e` the effort across the network and
# `power.f` the flow into it. MTKSL counts its through-variables positive into each
# component, so the adapter's own through-variable is negated. When the MTKSL side
# supplies power, the power along the bond is simply negative.
#
# The mapping is research 11's (`.scratch/mtk-upgrade/research/11-coupling-with-mtk-stdlib.md`):
# effort = across-variable in the impedance-analogy domains (electrical, hydraulic);
# flow = rate of the across-variable in the mobility-analogy ones (mechanical).

"""
    ElectricalAdapter(; name)

Bond to an MTKSL electrical network connected between pins `p` and `n`:
`e = p.v - n.v` (voltage across the network) and `f` = current into the network at `p`.
Its bond conventionally points into the adapter; a grounded network uses
`connect(adapter.n, ground.g)`.
"""
@component function ElectricalAdapter(; name)
    @named power = PowerPort()
    @named p = Pin()
    @named n = Pin()
    System([power.e ~ p.v - n.v, power.f ~ -p.i, 0 ~ p.i + n.i], t; name, systems = [power, p, n])
end

"""
    TranslationalAdapter(; name)

Bond to an MTKSL `Translational` (velocity-based) network at `flange`: `e` = force on the
network, `f = flange.v`. Its bond conventionally points into the adapter.
"""
@component function TranslationalAdapter(; name)
    @named power = PowerPort()
    @named flange = Translational.MechanicalPort()
    System([power.e ~ -flange.f, power.f ~ flange.v], t; name, systems = [power, flange])
end

"""
    TranslationalPositionAdapter(; name, s = nothing)

Bond to an MTKSL `TranslationalPosition` network at `flange`: `e` = force on the network,
`f = ds/dt`. `s` is the flange's initial position; leave it `nothing` when the MTKSL side
already sets it. Its bond conventionally points into the adapter.
"""
@component function TranslationalPositionAdapter(; name, s = nothing)
    @named power = PowerPort()
    @named flange = TranslationalPosition.Flange()
    @variables s(t) = s [guess = 0.0]
    System([s ~ flange.s, D(s) ~ power.f, power.e ~ -flange.f], t; name, systems = [power, flange])
end

"""
    RotationalAdapter(; name, phi = nothing)

Bond to an MTKSL `Rotational` network at `flange` (fixed support): `e` = torque on the
network, `f = dφ/dt`. `phi` is the flange's initial angle; leave it `nothing` when the
MTKSL side already sets it. Its bond conventionally points into the adapter.
"""
@component function RotationalAdapter(; name, phi = nothing)
    @named power = PowerPort()
    @named flange = Rotational.Flange()
    @variables phi(t) = phi [guess = 0.0]
    System([phi ~ flange.phi, D(phi) ~ power.f, power.e ~ -flange.tau], t; name, systems = [power, flange])
end

"""
    HydraulicAdapter(; name)

Bond to an MTKSL `IsothermalCompressible` network at `port`: `e = port.p` and `f` = volume
flow into the network at the fluid's reference density, `-port.dm / port.ρ`. Connect the
`HydraulicFluid` in its own statement (`connect(fluid, src.port)`), not in the same
`connect` as the adapter. Its bond conventionally points into the adapter.
"""
@component function HydraulicAdapter(; name)
    @named power = PowerPort()
    @named port = IsothermalCompressible.HydraulicPort()
    System([power.e ~ port.p, power.f ~ -port.dm / port.ρ], t; name, systems = [power, port])
end
