# Freely-propagating premixed laminar flame, mirroring Cantera's `adiabatic_flame`
# example.  Computes the laminar flame speed of a stoichiometric methane/air
# mixture on the GRI-Mech 3.0 mechanism.

using Cantera

gas = Solution("gri30.yaml")
set_TPX!(gas, 300.0, one_atm, "CH4:1, O2:2, N2:7.52")

flame = FreeFlame(gas; width=0.03)
solve!(flame; auto=true)          # staged energy-off -> energy-on -> refine

println("laminar flame speed Su = ", flame_speed(flame), " m/s")
println("adiabatic flame T      = ", maximum(flame_T(flame)), " K")
println("grid points            = ", length(grid(flame)))

# Peak OH mole fraction through the flame front
oh = flame_X(flame, "OH")
println("peak X(OH)             = ", maximum(oh))
