"""
High-Pressure transport using two models
========================================

Two high-pressure gas transport models are available in Cantera: the 'high-pressure'
and the 'high-pressure-chung' models. These models utilize critical fluid properties
and other fluid-specific parameters to calculate transport properties at high pressures.
These models are useful for fluids that are supercritical where ideal gas assumptions
do not yield accurate results.

Requires:  cantera >= 3.1.0

.. tags:: Python, transport, high-pressure
"""
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

# Plot the difference between the high-pressure viscosity model and
# what an ideal gas model would give for Methane.
gas = ct.Solution('methane_co2.yaml')

pressures = np.linspace(101325, 6e7, 100)
# Collect ideal viscosities
ideal_viscosities = []
for pressure in pressures:
    gas.TPX = 350, pressure, 'CH4:1.0'
    ideal_viscosities.append(gas.viscosity)

# Collect real viscosities using the high-pressure-chung
gas.transport_model = 'high-pressure-chung'
real_viscosities_1 = []
for pressure in pressures:
    gas.TPX = 350, pressure, 'CH4:1.0'
    real_viscosities_1.append(gas.viscosity)

# Collect real viscosities using the high-pressure model
gas.transport_model = 'high-pressure'
real_viscosities_2 = []
for pressure in pressures:
    gas.TPX = 350, pressure, 'CH4:1.0'
    real_viscosities_2.append(gas.viscosity)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(pressures, ideal_viscosities, label='Ideal Viscosity')
plt.plot(pressures, real_viscosities_1, label='High-Pressure Viscosity (Chung)')
plt.plot(pressures, real_viscosities_2, label='High-Pressure Viscosity (Lucas)')
plt.xlabel('Pressure (Pa)')
plt.ylabel('Viscosity  (Pa·s)')
plt.title('Methane Viscosity Model Comparison')
plt.legend()
plt.grid(True)
plt.show()



