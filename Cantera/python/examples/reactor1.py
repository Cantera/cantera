"""

  Constant-pressure, adiabatic kinetics simulation.

"""
from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from Cantera import rxnpath

gri3 = GRI30()

gri3.setState_TPX(1001.0, OneAtm, 'H2:2,O2:1,N2:4')
r   = Reactor(gri3)

env = Reservoir(Air())

# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = Wall(r,env)
w.set(K = 1.0e6)   # set expansion parameter. dV/dt = K(P_1 - P_2)
w.set(A = 1.0)

time = 0.0
for n in range(100):
    time += 1.e-5
    r.advance(time)
    env.advance(time)    
    print '%10.3e %10.3f %10.3f %14.6e' % (r.time(), r.temperature(), 
                                           r.pressure(), r.intEnergy_mass())

#print gri3

