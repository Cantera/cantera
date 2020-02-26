# coding: utf-8

"""
Simulate two counter-flow jets of reactants shooting into each other. This
simulation differs from the similar premixed_counterflow_flame.py example as the
latter simulates a jet of reactants shooting into products.

Requires: cantera >= 2.5.0
"""

import cantera as ct
import numpy as np
import sys


# Differentiation function for data that has variable grid spacing Used here to
# compute normal strain-rate
def derivative(x, y):
    dydx = np.zeros(y.shape, y.dtype.type)

    dx = np.diff(x)
    dy = np.diff(y)
    dydx[0:-1] = dy/dx

    dydx[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])

    return dydx


def computeStrainRates(oppFlame):
    # Compute the derivative of axial velocity to obtain normal strain rate
    strainRates = derivative(oppFlame.grid, oppFlame.velocity)

    # Obtain the location of the max. strain rate upstream of the pre-heat zone.
    # This is the characteristic strain rate
    maxStrLocation = abs(strainRates).argmax()
    minVelocityPoint = oppFlame.velocity[:maxStrLocation].argmin()

    # Characteristic Strain Rate = K
    strainRatePoint = abs(strainRates[:minVelocityPoint]).argmax()
    K = abs(strainRates[strainRatePoint])

    return strainRates, strainRatePoint, K


def computeConsumptionSpeed(oppFlame):

    Tb = max(oppFlame.T)
    Tu = min(oppFlame.T)
    rho_u = max(oppFlame.density)

    integrand = oppFlame.heat_release_rate/oppFlame.cp

    total_heat_release = np.trapz(integrand, oppFlame.grid)
    Sc = total_heat_release/(Tb - Tu)/rho_u

    return Sc


# This function is called to run the solver
def solveOpposedFlame(oppFlame, massFlux=0.12, loglevel=1,
                      ratio=2, slope=0.3, curve=0.3, prune=0.05):
    """
    Execute this function to run the Oppposed Flow Simulation This function
    takes a CounterFlowTwinPremixedFlame object as the first argument
    """

    oppFlame.reactants.mdot = massFlux
    oppFlame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

    oppFlame.show_solution()
    oppFlame.solve(loglevel, auto=True)

    # Compute the strain rate, just before the flame. This is not necessarily
    # the maximum We use the max. strain rate just upstream of the pre-heat zone
    # as this is the strain rate that computations comprare against, like when
    # plotting Su vs. K
    strainRates, strainRatePoint, K = computeStrainRates(oppFlame)

    return np.max(oppFlame.T), K, strainRatePoint


# Select the reaction mechanism
gas = ct.Solution('gri30.yaml')

# Create a CH4/Air premixed mixture with equivalence ratio=0.75, and at room
# temperature and pressure.
gas.set_equivalence_ratio(0.75, 'CH4', {'O2': 1.0, 'N2': 3.76})
gas.TP = 300, ct.one_atm

# Set the velocity
axial_velocity = 2.0  # in m/s

# Domain half-width of 2.5 cm, meaning the whole domain is 5 cm wide
width = 0.025

# Done with initial conditions
# Compute the mass flux, as this is what the Flame object requires
massFlux = gas.density * axial_velocity  # units kg/m2/s

# Create the flame object
oppFlame = ct.CounterflowTwinPremixedFlame(gas, width=width)

# Uncomment the following line to use a Multi-component formulation. Default is
# mixture-averaged
# oppFlame.transport_model = 'Multi'

# Now run the solver. The solver returns the peak temperature, strain rate and
# the point which we ascribe to the characteristic strain rate.

(T, K, strainRatePoint) = solveOpposedFlame(oppFlame, massFlux, loglevel=1)

# You can plot/see all state space variables by calling oppFlame.foo where foo
# is T, Y[i], etc. The spatial variable (distance in meters) is in oppFlame.grid
# Thus to plot temperature vs distance, use oppFlame.grid and oppFlame.T

Sc = computeConsumptionSpeed(oppFlame)

print("Peak temperature: {0:.1f} K".format(T))
print("Strain Rate: {0:.1f} 1/s".format(K))
print("Consumption Speed: {0:.2f} cm/s".format(Sc*100))
oppFlame.write_csv("premixed_twin_flame.csv", quiet=False)

# Generate plots to see results, if user desires
if '--plot' in sys.argv:

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6), facecolor='white')

    # Axial Velocity Plot
    plt.subplot(1, 2, 1)
    plt.plot(oppFlame.grid, oppFlame.velocity, 'r', lw=2)
    plt.xlim(oppFlame.grid[0], oppFlame.grid[-1])
    plt.xlabel('Distance (m)')
    plt.ylabel('Axial Velocity (m/s)')

    # Identify the point where the strain rate is calculated
    plt.plot(oppFlame.grid[strainRatePoint],
             oppFlame.velocity[strainRatePoint], 'gs')
    plt.annotate('Strain-Rate point',
                 xy=(oppFlame.grid[strainRatePoint],
                     oppFlame.velocity[strainRatePoint]),
                 xytext=(0.001, 0.1),
                 arrowprops={'arrowstyle': '->'})

    # Temperature Plot
    plt.subplot(1, 2, 2)
    plt.plot(oppFlame.grid, oppFlame.T, 'b', lw=2)
    plt.xlim(oppFlame.grid[0], oppFlame.grid[-1])
    plt.xlabel('Distance (m)')
    plt.ylabel('Temperature (K)')

    plt.tight_layout()
    plt.show()

else:
    print('************')
    print('Plotting option not enabled. Re-run script with --plot to see key plots.')
    print('************')
