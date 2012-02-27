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

   The two volumes are connected by an adiabatic free piston.
   The piston speed is proportional to the pressure difference
   between the two chambers.

   Note that each side uses a *different* reaction mechanism

"""
from Cantera import *
from Cantera.Reactor import *
import sys

fmt = '%10.3f  %10.1f  %10.4f  %10.4g  %10.4g  %10.4g  %10.4g'
print '%10s  %10s  %10s  %10s  %10s  %10s %10s' % ('time [s]','T1 [K]','T2 [K]',
                                              'V1 [m^3]', 'V2 [m^3]',
                                              'V1+V2 [m^3]','X(CO)')

gas1 = importPhase('h2o2.cti')
gas1.set(T = 900.0, P = OneAtm, X = 'H2:2, O2:1, AR:20')

gas2 = GRI30()
gas2.set(T = 900.0, P = OneAtm, X = 'CO:2, H2O:0.01, O2:5')

r1 = Reactor(gas1, volume = 0.5)
r2 = Reactor(gas2, volume = 0.1)
w = Wall(left = r1, right = r2, K = 1.0e3)

reactors = ReactorNet([r1, r2])

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
    reactors.advance(time)
    print fmt % (time, r1.temperature(), r2.temperature(),
                 r1.volume(), r2.volume(), r1.volume() + r2.volume(),
                 r2.moleFraction('CO'))

    tim.append(time)
    t1.append(r1.temperature())
    t2.append(r2.temperature())
    v1.append(r1.volume())
    v2.append(r2.volume())
    v.append(r1.volume() + r2.volume())
    xco.append(r2.moleFraction('CO'))
    xh2.append(r1.moleFraction('H2'))


# plot the results if matplotlib is installed.
# see http://matplotlib.sourceforge.net to get it
args = sys.argv
if len(args) > 1 and (args[1] == '-plot' or
                      args[1] == '-p' or
                      args[1] == '--plot'):
    try:
        from matplotlib.pylab import *
        clf
        subplot(2,2,1)
        plot(tim,t1,'-',tim,t2,'r-')
        xlabel('Time (s)');
        ylabel('Temperature (K)');
        subplot(2,2,2)
        plot(tim,v1,'-',tim,v2,'r-',tim,v,'g-')
        xlabel('Time (s)');
        ylabel('Volume (m3)');
        subplot(2,2,3)
        plot(tim,xco);
        xlabel('Time (s)');
        ylabel('CO Mole Fraction (right)');
        subplot(2,2,4)
        plot(tim,xh2);
        xlabel('Time (s)');
        ylabel('H2 Mole Fraction (left)');
        show()
    except:
        print """matplotlib required.
        http://matplotlib.sourceforge.net"""

else:
    print """To view a plot of these results, run this script with the option -plot"""
