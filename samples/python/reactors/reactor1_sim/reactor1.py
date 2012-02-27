"""

  Constant-pressure, adiabatic kinetics simulation.

"""
import sys

from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from Cantera import rxnpath

gri3 = GRI30()

gri3.set(T = 1001.0, P = OneAtm, X = 'H2:2,O2:1,N2:4')
r   = Reactor(gri3)

env = Reservoir(Air())

# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = Wall(r,env)
w.set(K = 1.0e6)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.set(A = 1.0)

sim = ReactorNet([r])
time = 0.0
tim = zeros(100,'d')
data = zeros([100,5],'d')

for n in range(100):
    time += 1.e-5
    sim.advance(time)
    tim[n] = time
    data[n,0] = r.temperature()
    data[n,1] = r.moleFraction('OH')
    data[n,2] = r.moleFraction('H')
    data[n,3] = r.moleFraction('H2')
    print '%10.3e %10.3f %10.3f %14.6e' % (sim.time(), r.temperature(),
                                           r.pressure(), r.intEnergy_mass())


# plot the results if matplotlib is installed.
# see http://matplotlib.sourceforge.net to get it
args = sys.argv
if len(args) > 1 and args[1] == '-plot':
    try:
        from matplotlib.pylab import *
        clf
        subplot(2,2,1)
        plot(tim,data[:,0])
        xlabel('Time (s)');
        ylabel('Temperature (K)');
        subplot(2,2,2)
        plot(tim,data[:,1])
        xlabel('Time (s)');
        ylabel('OH Mole Fraction');
        subplot(2,2,3)
        plot(tim,data[:,2]);
        xlabel('Time (s)');
        ylabel('H Mole Fraction');
        subplot(2,2,4)
        plot(tim,data[:,3]);
        xlabel('Time (s)');
        ylabel('H2 Mole Fraction');
        show()
    except:
        pass
else:
    print """To view a plot of these results, run this script with the option -plot"""
