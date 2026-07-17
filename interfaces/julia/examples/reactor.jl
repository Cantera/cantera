# Constant-volume adiabatic ignition of a hydrogen/air mixture, integrated with
# a ReactorNet.  Mirrors Cantera's canonical `reactor1.py` example.

using Cantera

gas = Solution("gri30.yaml")
set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")

reactor = IdealGasReactor(gas)
net = ReactorNet(reactor)

println("# t [s]    T [K]      p [Pa]")
t = 0.0
tend = 1e-3
while t < tend
    global t = step!(net)
    println(rpad(round(t, sigdigits=4), 10), " ",
            rpad(round(temperature(reactor), digits=2), 10), " ",
            round(pressure(reactor), digits=1))
end

println("\nFinal temperature: ", temperature(reactor), " K")
