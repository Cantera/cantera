"""
Constant-pressure, adiabatic kinetics simulation with sensitivity analysis
"""

import sys
import numpy as np

import cantera as ct

gri3 = ct.Solution('gri30.xml')
temp = 1500.0
pres = ct.one_atm

gri3.TPX = temp, pres, 'CH4:0.1, O2:2, N2:7.52'
r = ct.IdealGasConstPressureReactor(gri3, name='R1')
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


n_times = 400
tim = np.zeros(n_times)
data = np.zeros((n_times,6))

time = 0.0
for n in range(n_times):
    time += 5.0e-6
    sim.advance(time)
    tim[n] = 1000 * time
    data[n,0] = r.T
    data[n,1:4] = r.thermo['OH','H','CH4'].X

    # sensitivity of OH to reaction 2
    data[n,4] = sim.sensitivity('OH',2)

    # sensitivity of OH to reaction 3
    data[n,5] = sim.sensitivity('OH',3)

    print('%10.3e %10.3f %10.3f %14.6e %10.3f %10.3f' %
          (sim.time, r.T, r.thermo.P, r.thermo.u,  data[n,4],  data[n,5]))

# plot the results if matplotlib is installed.
# see http://matplotlib.org/ to get it
if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    plt.subplot(2,2,1)
    plt.plot(tim,data[:,0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2,2,2)
    plt.plot(tim,data[:,1])
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2,2,3)
    plt.plot(tim,data[:,2])
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2,2,4)
    plt.plot(tim,data[:,3])
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction')
    plt.tight_layout()

    plt.figure(2)
    plt.plot(tim,data[:,4],'-',tim,data[:,5],'-g')
    plt.legend([sim.sensitivity_parameter_name(2),sim.sensitivity_parameter_name(3)],'best')
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Sensitivity')
    plt.tight_layout()
    plt.show()
else:
    print("""To view a plot of these results, run this script with the option '--plot""")
