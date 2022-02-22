"""
This example demonstrates how to set a mixture according to equivalence ratio
and mixture fraction.

Requires: cantera >= 2.6.0
"""

import cantera as ct

gas = ct.Solution('gri30.yaml')

# Set the mixture composition according to the stoichiometric mixture
# (equivalence ratio phi = 1). The fuel composition in this example
# is set to 100 mol-% CH4 and the oxidizer to 21 mol-% O2 and
# 79 mol-% N2.
gas.set_equivalence_ratio(phi=1.0, fuel="CH4:1", oxidizer="O2:0.21,N2:0.79")

# The following is a short hand notation and equivalent to the call above
gas.set_equivalence_ratio(1.0, "CH4", "O2:0.21,N2:0.79")

# By default, the composition of fuel and oxidizer are interpreted as mole
# fractions. If the compositions are given in mass fractions, an
# additional argument can be provided. Here, the fuel is 100 mass-% CH4
# and the oxidizer is 23.3 mass-% O2 and 76.7 mass-% N2.
gas.set_equivalence_ratio(1.0, "CH4", "O2:0.233,N2:0.767", basis='mass')

# This function can be used to compute the equivalence ratio for any mixture.
# The first two arguments specify the compositions of the fuel and oxidizer.
# An optional third argument "basis" indicates if fuel and oxidizer compositions
# are provided in terms of mass or mole fractions. Default is mole fractions.
# Note that for all functions shown here, the compositions are normalized
# internally so the species fractions do not have to sum to unity
phi = gas.equivalence_ratio(fuel="CH4", oxidizer="O2:233,N2:767", basis='mass')
print("phi = {:1.3f}".format(phi))

# If the compositions of fuel and oxidizer are unknown, the function can
# be called without arguments. This assumes that all C, H and S atoms come from
# the fuel and all O atoms from the oxidizer. In this example, the fuel was set
# to be pre CH4 and the oxidizer O2:0.233,N2:0.767 so that the assumption is true
# and the same equivalence ratio as above is computed
phi = gas.equivalence_ratio()
print("phi = {:1.3f}".format(phi))

# Instead of working with equivalence ratio, mixture fraction can be used.
# The mixture fraction is always kg fuel / (kg fuel + kg oxidizer), independent
# of the basis argument. For example, the mixture fraction Z can be computed as
# follows. Again, the compositions by default are interpreted as mole fractions
Z = gas.mixture_fraction(fuel="CH4", oxidizer="O2:0.21,N2:0.79")
print("Z = {:1.3f}".format(Z))
# By default, the mixture fraction is the Bilger mixture fraction. Instead,
# a mixture fraction based on a single element can be used. In this example,
# the following two ways of computing Z are the same:
Z = gas.mixture_fraction(fuel="CH4", oxidizer="O2:0.21,N2:0.79", element="Bilger")
Z = gas.mixture_fraction(fuel="CH4", oxidizer="O2:0.21,N2:0.79", element="C")

# Since the fuel in this example is pure methane and the oxidizer is air,
# the mixture fraction is the same as the mass fraction of methane in the mixture
print("mass fraction of CH4 = {:1.3f}".format(gas["CH4"].Y[0]))

# To set a mixture according to the mixture fraction, the following function
# can be used. In this example, the final fuel/oxidizer mixture
# contains 5.5 mass-% CH4:
Z = gas.set_mixture_fraction(0.055, fuel="CH4", oxidizer="O2:0.21,N2:0.79")
print("mass fraction of CH4 = {:1.3f}".format(gas["CH4"].Y[0]))

# Mixture fraction and equivalence ratio are invariant to the reaction progress.
# For example, they stay constant if the mixture composition changes to the burnt
# state or for any intermediate state
fuel = "CH4:1"
oxidizer = "O2:0.21,N2:0.79"
gas.set_equivalence_ratio(1, fuel, oxidizer)
gas.equilibrate('HP')
phi_burnt = gas.equivalence_ratio(fuel, oxidizer)
Z_burnt = gas.mixture_fraction(fuel, oxidizer)
print("phi(burnt) = {:1.3f}".format(phi_burnt))
print("Z(burnt) = {:1.3f}".format(Z_burnt))

# If fuel and oxidizer compositions are specified consistently, then
# equivalence_ratio and set_equivalence_ratio are consistent as well, as
# shown in the following example
gas.set_equivalence_ratio(2.5,fuel="CH4:1,O2:0.01,CO:0.05,N2:0.1",
        oxidizer="O2:0.2,N2:0.8,CO2:0.05,CH4:0.01")

phi = gas.equivalence_ratio(fuel="CH4:1,O2:0.01,CO:0.05,N2:0.1",
        oxidizer="O2:0.2,N2:0.8,CO2:0.05,CH4:0.01")
print("phi = {:1.3f}".format(phi)) # prints 2.5

# Without specifying the fuel and oxidizer compositions, it is assumed that
# all C, H and S atoms come from the fuel and all O atoms from the oxidizer,
# which is not true for this example. Therefore, the following call gives a
# different equivalence ratio based on that assumption
phi = gas.equivalence_ratio()
print("phi = {:1.3f}".format(phi))

# After computing the mixture composition for a certain equivalence ratio given
# a fuel and mixture composition, the mixture can optionally be diluted. The
# following function will first create a mixture with equivalence ratio 2 from pure
# hydrogen and oxygen and then dilute it with H2O. In this example, the final mixture
# consists of 30 mol-% H2O and 70 mol-% of the H2/O2 mixture at phi=2
gas.set_equivalence_ratio(2.0, "H2:1", "O2:1", fraction={"diluent":0.3}, diluent="H2O")
print("mole fraction of H2O = {:1.3f}".format(gas["H2O"].X[0])) # final mixture contains 30 mol-% H2O
print("ratio of H2/O2: {:1.3f}".format(gas["H2"].X[0]/gas["O2"].X[0])) # according to phi=2

# Another option is to specify the fuel or oxidizer fraction in the final mixture.
# The following example creates a mixture with equivalence ratio 2 from pure
# hydrogen and oxygen (same as above) and then dilutes it with a mixture of 50 mass-%
# CO2 and 50 mass-% H2O so that the mass fraction of the fuel in the final mixture is 0.1
gas.set_equivalence_ratio(2.0, "H2", "O2", fraction={"fuel":0.1}, diluent="CO2:0.5,H2O:0.5", basis="mass")
print("mole fraction of H2O = {:1.3f}".format(gas["H2"].Y[0])) # final mixture contains 10 mass-% fuel

# To compute the equivalence ratio given a diluted mixture, a list of
# species names can be provided which will be considered for computing phi.
# In this example, the diluents H2O and CO2 are ignored and only H2 and O2 are
# considered to get the equivalence ratio
phi = gas.equivalence_ratio(fuel="H2", oxidizer="O2", include_species=["H2","O2"])
print("phi = {:1.3f}".format(phi)) # prints 2
