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
ar.setState_TPX(1000.0, 20.0*OneAtm, 'AR:1')

# create a reactor to represent the side of the cylinder filled with argon
r1   = Reactor(ar)


# create a reservoir for the environment, and fill it with air.
env = Reservoir(Air())


# use GRI-Mech 3.0 for the methane/air mixture, and set its initial state
gri3 = GRI30()

gri3.setState_TPX(500.0, 0.1*OneAtm, 'CH4:1.1, O2:2, N2:7.52')

# create a reactor for the methane/air side
r2 = Reactor(gri3)


#---------------------------------------------------------------------
# Now couple the reactors by defining common walls that may move (a piston)
# or conduct heat
#----------------------------------------------------------------------

# add a flexible wall (a piston) between r2 and r1
w = Wall(r2, r1)
w.set(area = 2.0, K=1.1e-4)


# heat loss to the environment. Heat loss always occur through walls,
# so we create a wall separating r1 from the environment, give it a
# non-zero area, and specify the overall heat transfer coefficient
# through the wall.
w2 = Wall(r1, env)
w2.set(area = 0.5, U=100.0)

# Now the problem is set up, and we're ready to solve it.
print 'finished setup, begin solution...'

time = 0.0
f = open('piston.csv','w')
writeCSV(f,['time (s)','T2 (K)','P2 (Pa)','V2 (m3)',
            'T1 (K)','P1 (Pa)','V1 (m3)'])
for n in range(300):
    time += 4.e-5
    print time, r2.temperature(),n
    r1.advance(time)    
    r2.advance(time)
    writeCSV(f, [r2.time(), r2.temperature(), r2.pressure(), r2.volume(),
                 r1.temperature(), r1.pressure(), r1.volume()])
f.close()
import os
print 'Output written to file piston.csv'
print 'Directory: '+os.getcwd()

