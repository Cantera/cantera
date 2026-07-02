# Reactor network with mass flow: a well-stirred-reactor-like setup where a hot
# reservoir feeds an IdealGasReactor through a MassFlowController and the reactor
# exhausts to a downstream reservoir through a Valve.  Demonstrates reservoirs,
# flow devices, and time integration.

using Cantera

# Upstream reservoir: fresh reactants.
upstream_gas = Solution("gri30.yaml")
set_TPX!(upstream_gas, 300.0, one_atm, "CH4:1, O2:2, N2:7.52")
upstream = Reservoir(upstream_gas)

# The reactor, initially filled with hot products to ignite the incoming flow.
reactor_gas = Solution("gri30.yaml")
set_TPX!(reactor_gas, 1800.0, one_atm, "CH4:1, O2:2, N2:7.52")
equilibrate!(reactor_gas, "TP")
reactor = IdealGasReactor(reactor_gas)

# Downstream exhaust reservoir.
downstream_gas = Solution("gri30.yaml")
set_TPX!(downstream_gas, 300.0, one_atm, "N2:1")
downstream = Reservoir(downstream_gas)

# Connect: fixed inlet mass flow, pressure-controlled outlet.
mfc = MassFlowController(upstream, reactor; mdot=0.05)
valve = Valve(reactor, downstream; K=1.0)

net = ReactorNet(reactor)
set_tolerances!(net; rtol=1e-8, atol=1e-14)

println("# t [s]      T [K]       inlet mdot [kg/s]")
for t in range(0, 2e-2; length=6)
    advance!(net, t)
    println(rpad(round(t, sigdigits=3), 10), "  ",
            rpad(round(temperature(reactor), digits=1), 10), "  ",
            round(mass_flow_rate(mfc), sigdigits=4))
end
