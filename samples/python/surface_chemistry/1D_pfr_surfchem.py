"""
A 1-D steady state plug-flow reactor with surface chemistry.

Assumes an adiabatic, frictionless, constant area, cylindrical reactor.

Based off the jupyter notebook created by Yuanjie Jiang.

Requires: cantera >= 3.0, matplotlib >= 2.0
Keywords: catalysis, plug flow reactor, reactor network, surface chemistry
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

mech = 'SiF4_NH3_mec.yaml'

# import the model for gas-Si-N interface and adjacent bulk phases
gas_si_n_interface = ct.Interface(mech, 'SI3N4')
gas = gas_si_n_interface.adjacent['gas']

# Set the initial conditions
T0 = 1750  # K
p0 = 5000 * ct.one_atm / 760.0  # Pa ~2Torr
gas.TPX = T0, p0, "NH3:1, SiF4:1, N2:4, H2:3"
gas_si_n_interface.TP = T0, p0
D = 5.08e-2  # diameter of the tube [m]
Ac = np.pi * D**2 / 4  # cross section of the tube [m]
u0 = 11.53  # m/s initial velocity of the flow

reactor = ct.FlowReactor(gas, area=Ac)
reactor.mass_flow_rate = gas.density * u0 * Ac

rsurf = ct.ReactorSurface(gas_si_n_interface, reactor)
net = ct.ReactorNet([reactor])
soln = ct.SolutionArray(gas, extra=['x', 'speed', 'surf_coverages'])

# Integrate the reactor network. Note that "time" is really the coordinate along
# the length of the plug flow reactor
while net.time < 1.2:
    print(net.time, rsurf.coverages)
    net.step()
    soln.append(TDY=reactor.thermo.TDY,
                x=net.time,
                speed=reactor.speed,
                surf_coverages=rsurf.coverages)

# Plot the results
plt.rcParams['figure.constrained_layout.use'] = True
f, ax = plt.subplots(2, 3, figsize=(9,6))

# plot the temperature profile along the flow direction
ax[0, 0].plot(soln.x, soln.T[:], color='C3')
ax[0, 0].set(xlabel='Distance (m)', ylabel='Temperature (K)', title='Bulk properties')

# plot the pressure of the gas along the flow direction
ax_p = ax[0,0].twinx()
ax_p.plot(soln.x, soln.P, color='C2')
ax_p.set(ylabel='Pressure (Pa)')
ax_p.ticklabel_format(axis='y', style='sci',
                          scilimits=(-2, 2))  # scientific notation

# plot gas velocity along the flow direction
h_vel = ax[1, 0].plot(soln.x, soln.speed, color='C0', label='velocity')
ax[1, 0].set(xlabel='Distance (m)', ylabel='Velocity (m/s)')

# plot gas density along the flow direction
ax_rho = ax[1,0].twinx()
h_rho = ax_rho.plot(soln.x, soln.density, color='C1', label='density')
ax_rho.set(ylabel=r'Density ($\mathregular{kg/m^3}$)')
ax_rho.legend(handles=h_vel+h_rho)

# plot major and minor gas species separately
minor_idx = []
major_idx = []
for i, name in enumerate(gas.species_names):
    mean = np.mean(soln(name).Y)
    if mean >= 0.001:
        major_idx.append(i)
    elif mean >= 1e-10:
        minor_idx.append(i)

# plot major gas species along the flow direction
for j in major_idx:
    ax[0, 1].plot(soln.x, soln.Y[:,j], label=gas.species_name(j))
ax[0, 1].legend(fontsize=8, loc='best')
ax[0, 1].set(xlabel='Distance (m)', ylabel='Mass Fraction',
             title='Gas phase major species')

# plot minor gas species along the flow direction
for k, i in enumerate(minor_idx):
    style = '-' if k < 10 else '--'
    ax[1, 1].plot(soln.x, soln.Y[:,i], label=gas.species_name(i), linestyle=style)
ax[1, 1].legend(fontsize=7.5, loc='best')
ax[1, 1].set(xlabel='Distance (m)', ylabel='Mass Fraction',
             title='Gas phase minor species')

# plot the site fraction of the surface species along the flow direction
for i, name in enumerate(gas_si_n_interface.species_names):
    ax[0, 2].plot(soln.x, soln.surf_coverages[:, i], label=name)
ax[0, 2].legend(fontsize=8)
ax[0, 2].set(xlabel='Distance (m)', ylabel='Site Fraction', title='Surface species')

plt.show()
