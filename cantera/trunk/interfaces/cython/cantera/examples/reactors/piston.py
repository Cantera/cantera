"""
Gas 1: a stoichiometric H2/O2/Ar mixture
Gas 2: a wet CO/O2 mixture

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
"""

import sys

import cantera as ct

fmt = '%10.3f  %10.1f  %10.4f  %10.4g  %10.4g  %10.4g  %10.4g'
print('%10s  %10s  %10s  %10s  %10s  %10s %10s' % ('time [s]','T1 [K]','T2 [K]',
                                              'V1 [m^3]', 'V2 [m^3]',
                                              'V1+V2 [m^3]','X(CO)'))

gas1 = ct.Solution('h2o2.cti')
gas1.TPX = 900.0, ct.one_atm, 'H2:2, O2:1, AR:20'

gas2 = ct.Solution('gri30.xml')
gas2.TPX = 900.0, ct.one_atm, 'CO:2, H2O:0.01, O2:5'

r1 = ct.IdealGasReactor(gas1)
r1.volume = 0.5
r2 = ct.IdealGasReactor(gas2)
r2.volume = 0.1
w = ct.Wall(r1, r2, K=1.0e3)

net = ct.ReactorNet([r1, r2])

tim = []
t1 = []
t2 = []
v1 = []
v2 = []
v = []
xco = []
xh2 = []

for n in range(30):
    time = (n+1)*0.002
    net.advance(time)
    print(fmt % (time, r1.T, r2.T, r1.volume, r2.volume,
                 r1.volume + r2.volume, r2.thermo['CO'].X[0]))

    tim.append(time * 1000)
    t1.append(r1.T)
    t2.append(r2.T)
    v1.append(r1.volume)
    v2.append(r2.volume)
    v.append(r1.volume + r2.volume)
    xco.append(r2.thermo['CO'].X[0])
    xh2.append(r1.thermo['H2'].X[0])


# plot the results if matplotlib is installed.
if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    plt.subplot(2,2,1)
    plt.plot(tim,t1,'-',tim,t2,'r-')
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2,2,2)
    plt.plot(tim,v1,'-',tim,v2,'r-',tim,v,'g-')
    plt.xlabel('Time (ms)')
    plt.ylabel('Volume (m3)')
    plt.subplot(2,2,3)
    plt.plot(tim,xco)
    plt.xlabel('Time (ms)')
    plt.ylabel('CO Mole Fraction (right)')
    plt.subplot(2,2,4)
    plt.plot(tim,xh2)
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction (left)')
    plt.tight_layout()
    plt.show()

else:
    print("""To view a plot of these results, run this script with the option --plot""")
