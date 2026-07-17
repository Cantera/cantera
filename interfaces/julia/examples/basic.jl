# Basic usage of the Cantera Julia interface.
#
# Run from the repository root with:
#   CANTERA_LIBRARY_PATH=/path/to/lib julia --project=interfaces/julia \
#     interfaces/julia/examples/basic.jl

using Cantera

gas = Solution("gri30.yaml")
set_TPY!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")

println("T      = ", temperature(gas), " K")
println("p      = ", pressure(gas), " Pa")
println("rho    = ", density(gas), " kg/m^3")
println("cp     = ", cp_mass(gas), " J/kg/K")
println("lambda = ", thermal_conductivity(gas), " W/m/K")
println("mu     = ", viscosity(gas), " Pa*s")

wdot = net_production_rates(gas)
ich4 = species_index(gas, "CH4")
println("wdot(CH4) = ", wdot[ich4], " kmol/m^3/s")
