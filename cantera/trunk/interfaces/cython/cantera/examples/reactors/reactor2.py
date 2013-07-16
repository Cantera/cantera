"""
This script simulates the following situation. A closed cylinder with volume 2
m^3 is divided into two equal parts by a massless piston that moves with speed
proportional to the pressure difference between the two sides.  It is
initially held in place in the middle. One side is filled with 1000 K argon at
20 atm, and the other with a combustible 500 K methane/air mixture at 0.1 atm
(phi = 1.1). At t = 0 the piston is released and begins to move due to the
large pressure difference, compressing and heating the methane/air mixture,
which eventually explodes. At the same time, the argon cools as it expands.
The piston is adiabatic, but some heat is lost through the outer cylinder
walls to the environment.

Note that this simulation, being zero-dimensional, takes no account of shock
wave propagation. It is somewhat artifical, but nevertheless instructive.
"""

import sys
import os
import csv
import numpy as np

import cantera as ct

#-----------------------------------------------------------------------
# First create each gas needed, and a reactor or reservoir for each one.
#-----------------------------------------------------------------------

# create an argon gas object and set its state
ar = ct.Solution('argon.xml')
ar.TP = 1000.0, 20.0 * ct.one_atm

# create a reactor to represent the side of the cylinder filled with argon
r1 = ct.IdealGasReactor(ar)

# create a reservoir for the environment, and fill it with air.
env = ct.Reservoir(ct.Solution('air.xml'))

# use GRI-Mech 3.0 for the methane/air mixture, and set its initial state
gri3 = ct.Solution('gri30.xml')
gri3.TPX = 500.0, 0.2 * ct.one_atm, 'CH4:1.1, O2:2, N2:7.52'

# create a reactor for the methane/air side
r2 = ct.IdealGasReactor(gri3)

#-----------------------------------------------------------------------------
# Now couple the reactors by defining common walls that may move (a piston) or
# conduct heat
#-----------------------------------------------------------------------------

# add a flexible wall (a piston) between r2 and r1
w = ct.Wall(r2, r1, A=1.0, K=0.5e-4, U=100.0)

# heat loss to the environment. Heat loss always occur through walls, so we
# create a wall separating r1 from the environment, give it a non-zero area,
# and specify the overall heat transfer coefficient through the wall.
w2 = ct.Wall(r2, env, A=1.0, U=500.0)

sim = ct.ReactorNet([r1, r2])

# Now the problem is set up, and we're ready to solve it.
print('finished setup, begin solution...')

time = 0.0
n_steps = 300
outfile = open('piston.csv', 'w')
csvfile = csv.writer(outfile)
csvfile.writerow(['time (s)','T1 (K)','P1 (Bar)','V1 (m3)',
                  'T2 (K)','P2 (Bar)','V2 (m3)'])
temp = np.zeros((n_steps, 2))
pres = np.zeros((n_steps, 2))
vol = np.zeros((n_steps, 2))
tm = np.zeros(n_steps)

for n in range(n_steps):
    time += 4.e-4
    print(n, time, r2.T)
    sim.advance(time)
    tm[n] = time
    temp[n,:] = r1.T, r2.T
    pres[n,:] = 1.0e-5*r1.thermo.P, 1.0e-5*r2.thermo.P
    vol[n,:] = r1.volume, r2.volume
    csvfile.writerow([tm[n], temp[n,0], pres[n,0], vol[n,0],
                      temp[n,1], pres[n,1], vol[n,1]])
outfile.close()
print('Output written to file piston.csv')
print('Directory: '+os.getcwd())

if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2,2,1)
    h = plt.plot(tm, temp[:,0],'g-',tm, temp[:,1],'b-')
    #plt.legend(['Reactor 1','Reactor 2'],2)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')

    plt.subplot(2,2,2)
    plt.plot(tm, pres[:,0],'g-',tm, pres[:,1],'b-')
    #plt.legend(['Reactor 1','Reactor 2'],2)
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Bar)')

    plt.subplot(2,2,3)
    plt.plot(tm, vol[:,0],'g-',tm, vol[:,1],'b-')
    #plt.legend(['Reactor 1','Reactor 2'],2)
    plt.xlabel('Time (s)')
    plt.ylabel('Volume (m$^3$)')

    plt.figlegend(h, ['Reactor 1', 'Reactor 2'], loc='lower right')
    plt.tight_layout()
    plt.show()
else:
    print("""To view a plot of these results, run this script with the option -plot""")
