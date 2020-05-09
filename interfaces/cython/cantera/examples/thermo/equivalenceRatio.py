"""
This example demonstrates how to set a mixture according to equivlance ratio
and mixture fraction.

Requires: cantera >= 2.5.0
"""

import cantera as ct

gas = ct.Solution('gri30.yaml')

# fuel and oxidizer compositions
fuel = "CH4"
oxidizer = "O2:0.21,N2:0.79"

gas.TP = 300, ct.one_atm

# set the mixture composition according to the stoichiometric mixture (equivlance ratio = 1)
gas.set_equivalence_ratio(1, fuel, oxidizer)

# this function can be used to compute the equivalence ratio for any mixture
# an optional argument "basis" indicates if fuel and oxidizer compositions are
# provided in terms of mass or mole fractions. Default is mole fractions. In
# fuel and oxidizer are given in mass fractions, use basis='mass'
phi = gas.get_equivalence_ratio(fuel, oxidizer)
print("phi", phi)

# The equivalence ratio can also be computed from the elemental composition
# assuming that there is no oxygen in the fuel and no C,H and S elements in the oxidizer
phi = gas.get_equivalence_ratio()

# in this example, the result is the same as above
print("phi", phi)

# the mixture fraction Z can be computed as follows:
Z = gas.get_mixture_fraction(fuel, oxidizer)
print("Z", Z)

# The mixture fraction is kg fuel / kg oxidizer. Since the fuel in this example
# is pure methane and the oxidizer is air, the mixture fraction is the same as
# the mass fraction of methane in the mixture
print("mass fraction of CH4: ", gas["CH4"].Y)

# Mixture fraction and equivalence ratio are invariant to the reaction progress.
# For example, they stay constant if the mixture composition changes to the burnt
# state
gas.equilibrate('HP')
phi_burnt = gas.get_equivalence_ratio()
Z_burnt = gas.get_mixture_fraction(fuel, oxidizer)
print("phi(burnt):", phi_burnt)
print("Z(burnt):", Z)
