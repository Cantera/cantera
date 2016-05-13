# coding: utf-8

"""
Simulate two counter-flow jets of reactants shooting into each other. This
simulation differs from the similar premixed_counterflow_flame.py example as the
latter simulates a jet of reactants shooting into products.
"""

import cantera as ct
import numpy as np

# This function is called to run the solver
def solveOpposedFlame(oppFlame, massFlux=0.12, loglevel=1,
                      ratio=3, slope=0.15, curve=0.25, prune=0.05):
    """
    Execute this function to run the Oppposed Flow Simulation This function
    takes a CounterFlowTwinPremixedFlame object as the first argument
    """

    oppFlame.reactants.mdot = massFlux
    oppFlame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

    oppFlame.show_solution()
    oppFlame.solve(loglevel, auto=True)

    # Compute the strain rate, just before the flame. It also turns out to the
    # maximum. This is the strain rate that computations comprare against, like
    # when plotting Su vs. K
    peakStrain = np.max(np.gradient(oppFlame.u, np.gradient(oppFlame.grid)))
    return np.max(oppFlame.T), peakStrain


# Use the standard GRI3.0 Mechanism for CH4
gas = ct.Solution('gri30.cti')

# Create a CH4/Air premixed mixture with equivalence ratio=0.75, and at room
# temperature and pressure.
gas.set_equivalence_ratio(0.75, 'CH4', {'O2':1.0, 'N2':3.76})
gas.TP = 300, ct.one_atm

# Set the velocity
axial_velocity = 0.25 # in m/s

# Domain half-width of 2.5 cm, meaning the whole domain is 5 cm wide
width = 0.025

# Done with initial conditions

# Compute the mass flux, as this is what the Flame object requires
massFlux = gas.density * axial_velocity # units kg/m2/s
# Create the flame object
oppFlame = ct.CounterflowTwinPremixedFlame(gas, width=width)

# Now run the solver

# The solver returns the peak temperature and strain rate. You can plot/see all
# state space variables by calling oppFlame.foo where foo is T, Y[i] or whatever
# The spatial variable (distance in meters) is in oppFlame.grid Thus to plot
# temperature vs distance, use oppFlame.grid and oppFlame.T
(T, K) = solveOpposedFlame(oppFlame, massFlux)

print("Peak temperature: {0}".format(T))
print("Strain Rate: {0}".format(K))
oppFlame.write_csv("premixed_twin_flame.csv", quiet=False)
