from Cantera import *
import math

g, dbulk = importPhases('diamond.cti',['gas','diamond'])

d = importInterface('diamond.cti','diamond_100',phases = [g, dbulk])
ns = d.nSpecies()
mw = dbulk.molarMasses()[0]

t = 1200.0
x = g.moleFractions()
p = 20.0*OneAtm/760.0
g.setState_TPX(t, p, x)
ih = g.speciesIndex('H')

xh0 = x[ih]
f = open('d.csv','w')
#writeCSV(f, ['H mole Fraction', 'Rate']+d.speciesNames())
for n in range(20):
    x[ih] /= 1.4
    g.setState_TPX(t, p, x)
    d.advanceCoverages(10.0)
    cdot = d.netProductionRates(phase = dbulk)[0]
    mdot = mw*cdot
    rate = mdot/dbulk.density()
    print x[ih], rate*1.0e6*3600.0
    writeCSV(f,[x[ih],rate]+list(d.coverages()))
f.close()
