"""
Reactors separated by a moving piston
=====================================

Two reactors separated by a piston that moves with a speed proportional to the pressure
difference between the reactors.

- Gas 1: a stoichiometric H2/O2/Ar mixture
- Gas 2: a wet CO/O2 mixture

.. code:: none

    -------------------------------------
    |          ||                       |
    |          ||                       |
    |  gas 1   ||        gas 2          |
    |          ||                       |
    |          ||                       |
    -------------------------------------


The two volumes are connected by an adiabatic free piston. The piston speed is
proportional to the pressure difference between the two chambers.

Note that each side uses a *different* reaction mechanism

Requires: cantera >= 3.2.0, matplotlib >= 2.0

.. tags:: Python, combustion, reactor network, plotting
"""

import cantera as ct
import matplotlib.pyplot as plt
plt.rcParams['figure.constrained_layout.use'] = True

# %%
# Create objects representing the gases and reactors

gas1 = ct.Solution('h2o2.yaml')
gas1.TPX = 900.0, ct.one_atm, 'H2:2, O2:1, AR:20'

gas2 = ct.Solution('gri30.yaml')
gas2.TPX = 900.0, ct.one_atm, 'CO:2, H2O:0.01, O2:5'

r1 = ct.IdealGasReactor(gas1, clone=True)
r1.volume = 0.5
r2 = ct.IdealGasReactor(gas2, clone=True)
r2.volume = 0.1

# %%
# The wall is held fixed until t = 0.1 s, then released to allow the pressure to
# equilibrate.
def v(t):
    if t < 0.1:
        return 0.0
    else:
        return (r1.phase.P - r2.phase.P) * 1e-4

w = ct.Wall(r1, r2, velocity=v)

net = ct.ReactorNet([r1, r2])

# %%
# Run the simulation and collect the states of each reactor
states1 = ct.SolutionArray(r1.phase, extra=['t', 'volume'])
states2 = ct.SolutionArray(r2.phase, extra=['t', 'volume'])

fmt = '{:10.3f}  {:10.1f}  {:10.4f}  {:10.4g}  {:10.4g}  {:10.4g}  {:10.4g}'
print('{:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}'.format(
    'time [s]', 'T1 [K]', 'T2 [K]', 'V1 [m^3]', 'V2 [m^3]', 'Vtot [m^3]', 'X(CO)'))
for n in range(200):
    time = (n+1)*0.001
    net.advance(time)
    if n % 4 == 3:
        print(fmt.format(time, r1.T, r2.T, r1.volume, r2.volume,
                         r1.volume + r2.volume, r2.phase['CO'].X[0]))

    states1.append(r1.phase.state, t=1000*time, volume=r1.volume)
    states2.append(r2.phase.state, t=1000*time, volume=r2.volume)

# %%
# Plot the results
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(states1.t, states1.T, '-', states2.t, states2.T, 'r-')
ax[0,0].set(xlabel='Time (ms)', ylabel='Temperature (K)')

ax[0,1].plot(states1.t, states1.volume, '-', states2.t, states2.volume, 'r-',
             states1.t, states1.volume + states2.volume, 'g-')
ax[0,1].set(xlabel='Time (ms)', ylabel='Volume (m3)')

ax[1,0].plot(states2.t, states2('CO').X)
ax[1,0].set(xlabel='Time (ms)', ylabel='CO Mole Fraction (right)')

ax[1,1].plot(states1.t, states1('H2').X)
ax[1,1].set(xlabel='Time (ms)', ylabel='H2 Mole Fraction (left)')
plt.show()
