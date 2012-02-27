"""

  Constant-pressure, adiabatic kinetics simulation with sensitivity analysis

"""
import sys

from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *

gri3 = GRI30()
temp = 1500.0
pres = OneAtm

gri3.set(T = temp, P = pres, X = 'CH4:0.1, O2:2, N2:7.52')
r   = Reactor(gri3)

air = Air()
air.set(T = temp, P = pres)
env = Reservoir(air)

# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = Wall(r,env)
w.set(K = 1.0e6)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.set(A = 1.0)

# enable sensitivity with respect to the rates of the first 10
# reactions (reactions 0 through 9)
r.addSensitivityReaction(reactions = range(10))

sim = ReactorNet([r])

# set the tolerances for the solution and for the sensitivity
# coefficients
sim.setTolerances(rtol = 1.0e-6, atol = 1.0e-15,
                  rtolsens = 1.0e-5, atolsens = 1.0e-5)
time = 0.0
np = 400
tim = zeros(np,'d')
data = zeros([np,6],'d')

for n in range(np):
    time += 5.0e-6
    sim.advance(time)
    tim[n] = time
    data[n,0] = r.temperature()
    data[n,1] = r.moleFraction('OH')
    data[n,2] = r.moleFraction('H')
    data[n,3] = r.moleFraction('CH4')

    # sensitivity of OH to reaction 2
    data[n,4] = sim.sensitivity('OH',2)

    # sensitivity of OH to reaction 3
    data[n,5] = sim.sensitivity('OH',3)

    print '%10.3e %10.3f %10.3f %14.6e %10.3f %10.3f' % (sim.time(), r.temperature(),
                                           r.pressure(), r.intEnergy_mass(),  data[n,4],  data[n,5])


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
        figure(2)
        plot(tim,data[:,4],'-',tim,data[:,5],'-g')
        legend([r.sensParamName(2),r.sensParamName(3)],'best')
        xlabel('Time (s)');
        ylabel('OH Sensitivity');
        show()
    except:
        print 'could not make plots'
else:
    print """To view a plot of these results, run this script with the option -plot"""
