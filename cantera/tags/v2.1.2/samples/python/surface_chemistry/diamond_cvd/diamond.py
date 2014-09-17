# A CVD example. This example computes the growth rate of a diamond
# film according to a simplified version of a particular published
# growth mechanism (see file diamond.cti for details). Only the
# surface coverage equations are solved here; the gas composition is
# fixed. (For an example of coupled gas-phase and surface, see
# catcomb.py.)  Atomic hydrogen plays an important role in diamond
# CVD, and this example computes the growth rate and surface coverages
# as a function of [H] at the surface for fixed temperature and [CH3].

from Cantera import *
import math

print '\n\b******  CVD Diamond Example  ******\n'

# import the models for the gas and bulk diamond
g, dbulk = importPhases('diamond.cti',['gas','diamond'])

# import the model for the diamond (100) surface
d = importInterface('diamond.cti','diamond_100',phases = [g, dbulk])

ns = d.nSpecies()
mw = dbulk.molarMasses()[0]

t = 1200.0
x = g.moleFractions()
p = 20.0*OneAtm/760.0  # 20 Torr
g.set(T = t, P = p, X = x)

ih = g.speciesIndex('H')

xh0 = x[ih]
f = open('diamond.csv','w')
writeCSV(f, ['H mole Fraction', 'Growth Rate (microns/hour)']+d.speciesNames())
for n in range(20):
    x[ih] /= 1.4
    g.setState_TPX(t, p, x)
    d.advanceCoverages(10.0) # iintegrate the coverages to steady state
    carbon_dot = d.netProductionRates(phase = dbulk)[0]
    mdot = mw*carbon_dot
    rate = mdot/dbulk.density()
    writeCSV(f,[x[ih],rate*1.0e6*3600.0]+list(d.coverages()))
f.close()

print 'H concentration, growth rate, and surface coverages written to file diamond.csv'
