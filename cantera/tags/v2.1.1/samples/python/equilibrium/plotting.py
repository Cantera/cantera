
from matplotlib.pylab import *
from matplotlib import get_backend, interactive
import sys

# Print a warning if 'python' rather than 'pythonw' is invoked on a
# Mac.  Does nothing on other platforms.
def warnMac():
    if sys.platform == 'darwin':
        if sys.executable == '/usr/bin/python':
            b = get_backend()
            if b[-3:] == 'Agg':
                print 'Error: on a Mac, this script must be run with pythonw instead of python'
                print 'to display plots'
                return -1
    return 0


def plotEquilData(mix, phi, tad, xeq):

    if warnMac() < 0: return
    npoints = len(phi)

    nsp = mix.nSpecies()

    #titles = ['Major Species', 'Minor Species', 'N Minor Species']
    # assign species to 3 plots
    p = {}
    mm = 0
    for k in range(nsp):
        if amax(xeq[k,:]) > 0.01:
            p[k] = 0  # major species plot
        else:
            mm += 1
            if mix.nAtoms(k,'N') <= 0:
                p[k] = 1  # non-N minor species plot
            else:
                p[k] = 2  # N-containing minor species plot

    clf
    subplot(2,2,1)
    plot(phi,tad)
    xlabel('Equivalence Ratio')
    ylabel('Temperature (K)')
    title('Adiabatic Flame Temperature')
    axis([phi[1], phi[-1], amin(tad)-200.0, amax(tad)+200.0])

    # do three species plots
    for m in range(3):

        subplot(2,2,2+m);
        hold(True);

        for i in range(nsp):
            if p[i] == m:
                for j in range(npoints):
                    if xeq[i,j] <= 0.0:
                        xeq[i,j] = 1.0e-20
                if m == 0:
                    plot(phi, xeq[i,:])
                else:
                    semilogy(phi,xeq[i,:])

                xmax = amax(xeq[i,:])
                for j in range(npoints):
                    if xeq[i,j] == xmax: break
                offset = 0.0
                if j == 0:
                    offset = 0.3
                elif j >= npoints-1:
                    offset = -0.3
                text(phi[j]+offset,xeq[i,j],mix.speciesName(i))

        if m == 0:
            axis([phi[1], phi[-1], 0.0, 1.0]);
        else:
            axis([phi[1], phi[-1], 1.0e-14, 1]);

        xlabel('Equivalence Ratio');
        ylabel('Mole Fraction');
        #title(titles[m]);
        hold(False)

    show()
