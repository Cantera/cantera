"""
This example demonstrates how to set a mixture according to equivalence ratio
and mixture fraction.

Requires: cantera >= 2.5.0
"""

import cantera as ct

gas = ct.Solution('gri30.yaml')

# fuel and oxidizer compositions
fuel = "CH4"
oxidizer = "O2:0.21,N2:0.79"

gas.TP = 300, ct.one_atm

# set the mixture composition according to the stoichiometric mixture
# (equivalence ratio = 1)
gas.set_equivalence_ratio(1, fuel, oxidizer)

# This function can be used to compute the equivalence ratio for any mixture.
# An optional argument "basis" indicates if fuel and oxidizer compositions are
# provided in terms of mass or mole fractions. Default is mole fractions.
# If fuel and oxidizer are given in mass fractions, use basis='mass'
phi = gas.equivalence_ratio(fuel, oxidizer)
print("phi = {:1.3f}".format(phi))

# The equivalence ratio can also be computed from the elemental composition
# assuming that there is no oxygen in the fuel and no C,H and S elements
# in the oxidizer so that the composition of fuel and oxidizer can be omitted
phi = gas.equivalence_ratio()

# In this example, the result is the same as above
print("phi = {:1.3f}".format(phi))

# the mixture fraction Z can be computed as follows:
Z = gas.mixture_fraction(fuel, oxidizer)
print("Z = {:1.3f}".format(Z))

# The mixture fraction is kg fuel / (kg fuel + kg oxidizer). Since the fuel in
# this example is pure methane and the oxidizer is air, the mixture fraction
# is the same as the mass fraction of methane in the mixture
print("mass fraction of CH4 = {:1.3f}".format(gas["CH4"].Y[0]))

# Mixture fraction and equivalence ratio are invariant to the reaction progress.
# For example, they stay constant if the mixture composition changes to the burnt
# state
gas.equilibrate('HP')
phi_burnt = gas.equivalence_ratio(fuel, oxidizer)
Z_burnt = gas.mixture_fraction(fuel, oxidizer)
print("phi(burnt) = {:1.3f}".format(phi_burnt))
print("Z(burnt) = {:1.3f}".format(Z))
