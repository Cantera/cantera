"""
Constant-pressure, adiabatic kinetics simulation with sensitivity analysis

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import sys
import numpy as np

import cantera as ct

gas = ct.Solution('gri30.yaml')
temp = 1500.0
pres = ct.one_atm

gas.TPX = temp, pres, 'CH4:0.1, O2:2, N2:7.52'
r = ct.IdealGasConstPressureReactor(gas, name='R1')
sim = ct.ReactorNet([r])

# enable sensitivity with respect to the rates of the first 10
# reactions (reactions 0 through 9)
for i in range(10):
    r.add_sensitivity_reaction(i)

# set the tolerances for the solution and for the sensitivity coefficients
sim.rtol = 1.0e-6
sim.atol = 1.0e-15
sim.rtol_sensitivity = 1.0e-6
sim.atol_sensitivity = 1.0e-6

states = ct.SolutionArray(gas, extra=['t', 's2', 's3'])

for t in np.arange(0, 2e-3, 5e-6):
    sim.advance(t)
    s2 = sim.sensitivity('OH', 2)  # sensitivity of OH to reaction 2
    s3 = sim.sensitivity('OH', 3)  # sensitivity of OH to reaction 3
    states.append(r.thermo.state, t=1000*t, s2=s2, s3=s3)

    print('{:10.3e} {:10.3f} {:10.3f} {:14.6e} {:10.3f} {:10.3f}'.format(
        sim.time, r.T, r.thermo.P, r.thermo.u, s2, s3))

# plot the results if matplotlib is installed.
# see http://matplotlib.org/ to get it
if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    plt.subplot(2, 2, 1)
    plt.plot(states.t, states.T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(states.t, states('OH').X)
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(states.t, states('H').X)
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(states.t, states('CH4').X)
    plt.xlabel('Time (ms)')
    plt.ylabel('CH4 Mole Fraction')
    plt.tight_layout()

    plt.figure(2)
    plt.plot(states.t, states.s2, '-', states.t, states.s3, '-g')
    plt.legend([sim.sensitivity_parameter_name(2), sim.sensitivity_parameter_name(3)],
               'best')
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Sensitivity')
    plt.tight_layout()
    plt.show()
else:
    print("""To view a plot of these results, run this script with the option '--plot""")
