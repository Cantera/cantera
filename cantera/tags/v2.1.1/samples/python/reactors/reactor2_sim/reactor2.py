"""

This script simulates the following situation. A closed cylinder with
volume 2 m^3 is divided into two equal parts by a massless piston that
moves with speed proportional to the pressure difference between the
two sides.  It is initially held in place in the middle. One side is
filled with 1000 K argon at 20 atm, and the other with a combustible
500 K methane/air mixture at 0.1 atm (phi = 1.1). At t = 0 the piston
is released and begins to move due to the large pressure difference,
compressing and heating the methane/air mixture, which eventually
explodes. At the same time, the argon cools as it expands. The piston
is adiabatic, but some heat is lost through the outer cylinder walls
to the environment.

Note that this simulation, being zero-dimensional, takes no account of
shock wave propagation. It is somewhat artifical, but nevertheless
instructive.

"""

import sys
from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *

#-----------------------------------------------------------------------
# First create each gas needed, and a reactor or reservoir for each one.
#-----------------------------------------------------------------------

# create an argon gas object and set its state. This function is
# defined in module Cantera.gases, as are functions 'Air()', and
# 'GRI30()'

ar = Argon()
ar.set(T = 1000.0, P = 20.0*OneAtm, X = 'AR:1')

# create a reactor to represent the side of the cylinder filled with argon
r1   = Reactor(ar)


# create a reservoir for the environment, and fill it with air.
env = Reservoir(Air())


# use GRI-Mech 3.0 for the methane/air mixture, and set its initial state
gri3 = GRI30()

gri3.set(T = 500.0, P = 0.2*OneAtm, X = 'CH4:1.1, O2:2, N2:7.52')

# create a reactor for the methane/air side
r2 = Reactor(gri3)


#---------------------------------------------------------------------
# Now couple the reactors by defining common walls that may move (a piston)
# or conduct heat
#----------------------------------------------------------------------

# add a flexible wall (a piston) between r2 and r1
w = Wall(r2, r1)
w.set(area = 1.0, K=0.5e-4, U = 100.0)


# heat loss to the environment. Heat loss always occur through walls,
# so we create a wall separating r1 from the environment, give it a
# non-zero area, and specify the overall heat transfer coefficient
# through the wall.
w2 = Wall(r2, env)
w2.set(area = 1.0, U=500.0)

sim = ReactorNet([r1, r2])

# Now the problem is set up, and we're ready to solve it.
print 'finished setup, begin solution...'

time = 0.0
f = open('piston.csv','w')
writeCSV(f,['time (s)','T1 (K)','P1 (Bar)','V1 (m3)',
            'T2 (K)','P2 (Bar)','V2 (m3)'])
temp = zeros([300, 2], 'd')
pres = zeros([300, 2], 'd')
vol = zeros([300, 2], 'd')
tm = zeros(300,'d')
for n in range(300):
    time += 4.e-4
    print time, r2.temperature(),n
    sim.advance(time)
    tm[n] = time
    temp[n,:] = [r1.temperature(), r2.temperature()]
    pres[n,:] = [1.0e-5*r1.pressure(), 1.0e-5*r2.pressure()]
    vol[n,:] = [r1.volume(), r2.volume()]
    writeCSV(f, [tm[n], temp[n,0], pres[n,0], vol[n,0],
                 temp[n,1], pres[n,1], vol[n,1]])
f.close()
import os
print 'Output written to file piston.csv'
print 'Directory: '+os.getcwd()

args = sys.argv
if len(args) > 1 and args[1] == '-plot':
    try:
        from matplotlib.pylab import *
        clf
        subplot(2,2,1)
        plot(tm, temp[:,0],'g-',tm, temp[:,1],'b-')
        legend(['Reactor 1','Reactor 2'],2)
        xlabel('Time (s)');
        ylabel('Temperature (K)');

        subplot(2,2,2)
        plot(tm, pres[:,0],'g-',tm, pres[:,1],'b-')
        legend(['Reactor 1','Reactor 2'],2)
        xlabel('Time (s)');
        ylabel('Pressure (Bar)');

        subplot(2,2,3)
        plot(tm, vol[:,0],'g-',tm, vol[:,1],'b-')
        legend(['Reactor 1','Reactor 2'],2)
        xlabel('Time (s)');
        ylabel('Volume (m^3)');

        show()
    except:
        pass
else:
    print """To view a plot of these results, run this script with the option -plot"""
