---
file_format: mystnb
kernelspec:
  name: python3
---

```{py:currentmodule} cantera
```

# Calculating Heating Value of a Fuel

This guide demonstrates how to use Cantera's thermodynamic data to calculate the lower
[heating value](https://en.wikipedia.org/wiki/Heat_of_combustion) (LHV) and higher
heating value (HHV) of methane and other fuels.

## Heating value of Methane

The complete reaction for heating methane is:

$$\t{CH_4 + 2 O_2 \rightarrow CO_2 + 2 H_2O}$$

We compute the lower heating value (LHV) as the difference in enthalpy (per kg
*mixture*) between reactants and products at constant temperature and pressure, divided
by the mass fraction of fuel in the reactants.

```{code-cell} python
import cantera as ct
gas = ct.Solution("gri30.yaml")

# Set reactants state
gas.TPX = 298, 101325, "CH4:1, O2:2"
h1 = gas.enthalpy_mass
Y_CH4 = gas["CH4"].Y[0]  # returns an array, of which we only want the first element

# set state to complete combustion products without changing T or P
gas.TPX = None, None, "CO2:1, H2O:2"
h2 = gas.enthalpy_mass

LHV = -(h2 - h1) / Y_CH4 / 1e6
print(f"LHV = {LHV:.3f} MJ/kg")
```

The LHV is calculated assuming that water remains in the gas phase. However, more energy
can be extracted from the mixture if this water is condensed. This value is the higher
heating value (HHV).

The ideal gas mixture model used here cannot calculate this contribution directly.
However, Cantera also has a non-ideal equation of state which can be used to compute
this contribution.

```{code-cell} python
water = ct.Water()
# Set liquid water state, with vapor fraction x = 0
water.TQ = 298, 0
h_liquid = water.h
# Set gaseous water state, with vapor fraction x = 1
water.TQ = 298, 1
h_gas = water.h

# Calculate higher heating value
Y_H2O = gas["H2O"].Y[0]
HHV = -(h2 - h1 + (h_liquid - h_gas) * Y_H2O) / Y_CH4 / 1e6
print(f"HHV = {HHV:.3f} MJ/kg")
```

## Generalizing to arbitrary species

We can generalize this calculation by determining the composition of the products
automatically rather than directly specifying the product composition. This can be done
by computing the *elemental mole fractions* of the reactants mixture and noting that for
complete combustion, all of the carbon ends up as $\t{CO_2}$, all of the hydrogen ends
up as $\t{H_2O}$, and all of the nitrogen ends up as $\t{N_2}$. From this, we can
compute the ratio of these species in the products.

```{code-cell} python
def heating_value(fuel):
    """Returns the LHV and HHV for the specified fuel"""
    gas.TP = 298, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, "O2:1.0")
    h1 = gas.enthalpy_mass
    Y_fuel = gas[fuel].Y[0]

    # complete combustion products
    X_products = {
        "CO2": gas.elemental_mole_fraction("C"),
        "H2O": 0.5 * gas.elemental_mole_fraction("H"),
        "N2": 0.5 * gas.elemental_mole_fraction("N"),
    }

    gas.TPX = None, None, X_products
    Y_H2O = gas["H2O"].Y[0]
    h2 = gas.enthalpy_mass
    LHV = -(h2 - h1) / Y_fuel / 1e6
    HHV = -(h2 - h1 + (h_liquid - h_gas) * Y_H2O) / Y_fuel / 1e6
    return LHV, HHV


fuels = ["H2", "CH4", "C2H6", "C3H8", "NH3", "CH3OH"]
print("fuel   LHV (MJ/kg)   HHV (MJ/kg)")
for fuel in fuels:
    LHV, HHV = heating_value(fuel)
    print(f"{fuel:8s}   {LHV:7.3f}       {HHV:7.3f}")
```
