"""
Solution of a 1-D steady state plug-flow reactor with surface chemistry.
Assumes:
    - adiabatic
    - frictionless
    - constant area, cylindrical reactor

Based off the jupyter notebook created by Yuanjie Jiang
"""


import numpy as np
import cantera as ct
import matplotlib.pyplot as plt


mech = 'SiF4_NH3_mec.cti'
# import the models for gas and bulk
gas, bulk_si, bulk_n = ct.import_phases(mech, ['gas', 'SiBulk', 'NBulk'])
# import the model for gas-Si-N interface
gas_si_n_interface = ct.Interface(mech, 'SI3N4', [gas, bulk_si, bulk_n])

# Set the initial conditions
T0 = 1713  # K
p0 = 2 * ct.one_atm / 760.0  # Pa ~2Torr
gas.TPX = T0, p0, "NH3:6, SiF4:1"
bulk_si.TP = T0, p0
bulk_n.TP = T0, p0
gas_si_n_interface.TP = T0, p0
D = 5.08e-2  # diameter of the tube [m]
Ac = np.pi * D**2 / 4  # cross section of the tube [m]
mu = 5.7e-5  # kg/(m-s) dynamic viscosity
perim = np.pi * D  # perimeter of the tube
# calculate the site fractions of surface species at the entrance of the tube
# at steady state
u0 = 11.53  # m/s initial velocity of the flow


reactor = ct.FlowReactor(gas, area=Ac)
# set the mass flow-rate
reactor.mass_flow_rate = gas.density * u0 * Ac

rsurf = ct.ReactorSurface(gas_si_n_interface, reactor)
net = ct.ReactorNet([reactor])
soln = ct.SolutionArray(gas, extra=['t', 'speed', 'coverages'])

while net.time < 0.7:
    net.step()
    soln.append(T=reactor.thermo.T,
                P=reactor.thermo.P,
                X=reactor.thermo.X,
                t=net.time,
                speed=reactor.speed,
                coverages=rsurf.coverages)

f, ax = plt.subplots(4, 2, figsize=(9, 16), dpi=96)
# plot gas velocity along the flow direction
ax[0, 0].plot(soln.t, soln.speed, color='C0')
ax[0, 0].set_xlabel('Distance (m)')
ax[0, 0].set_ylabel('Velocity (m/s)')

# plot gas density along the flow direction
ax[0, 1].plot(soln.t, soln.density, color='C1')
ax[0, 1].set_xlabel('Distance (m)')
ax[0, 1].set_ylabel(r'Density ($\mathregular{kg/m^3}$)')
ax[0, 1].ticklabel_format(axis='y', style='sci',
                          scilimits=(-2, 2))  # scientific notation

# plot major and minor gas species separately
minor_idx = []
major_idx = []
for i, name in enumerate(gas.species_names):
    mean = np.mean(soln(name).Y)
    if mean <= 0.01:
        minor_idx.append(i)
    else:
        major_idx.append(i)

# plot minor gas species along the flow direction
for i in minor_idx:
    style = '-' if i < 10 else '--'
    ax[1, 0].plot(soln.t, soln(gas.species_names[i]).Y,
                  label=gas.species_names[i], linestyle=style)
ax[1, 0].legend(fontsize=7.5, loc='best')
ax[1, 0].set_xlabel('Distance (m)')
ax[1, 0].set_ylabel('Mass Fraction')
ax[1, 0].ticklabel_format(axis='y', style='sci',
                          scilimits=(-2, 2))  # scientific notation

# plot major gas species along the flow direction
for j in major_idx:
    ax[1, 1].plot(soln.t, soln(gas.species_names[j]).Y,
                  label=gas.species_names[j])
ax[1, 1].legend(fontsize=8, loc='best')
ax[1, 1].set_xlabel('Distance (m)')
ax[1, 1].set_ylabel('Mass Fraction')

# plot the pressure of the gas along the flow direction
ax[2, 0].plot(soln.t, soln.P, color='C2')
ax[2, 0].set_xlabel('Distance (m)')
ax[2, 0].set_ylabel('Pressure (Pa)')

# plot the site fraction of the surface species along the flow direction
for i, name in enumerate(gas_si_n_interface.species_names):
    ax[2, 1].plot(soln.t, soln.coverages[:, i], label=name)
ax[2, 1].legend(fontsize=8)
ax[2, 1].set_xlabel('Distance (m)')
ax[2, 1].set_ylabel('Site Fraction')

# plot the temperature profile along the flow direction
ax[3, 0].plot(soln.t, soln.T[:], color='C3')
ax[3, 0].set_xlabel('Distance (m)')
ax[3, 0].set_ylabel('Temperature (K)')
f.tight_layout(pad=0.5)
plt.show()
