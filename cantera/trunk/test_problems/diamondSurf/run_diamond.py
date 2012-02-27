from Cantera import *
# import the bulk phases
g, dbulk = importPhases('diamond.cti', ['gas','diamond'])

# import the interface
d = importInterface('diamond.cti', 'diamond_100', phases = [g, dbulk])

mw = dbulk.molarMasses()[0] #mol. wt. of carbin

t = g.temperature()
p = g.pressure()
x = g.moleFractions()
ih = g.speciesIndex('H')

f = open('d.csv', 'w')
for n in range (20):
    x[ih] /= 1.4
    g.setState_TPX(t, p, x)
    # integrate the coverage equations to steady state
    d.advanceCoverages(100.0)
    cdot = d.netProductionRates(phase = dbulk) [0] #net rate of C(d production /m^2
    mdot = mw * cdot
    linear_rate = mdot/dbulk.density()
    writeCSV(f, [x[ih], rate]+list(d.coverages()))
f.close()
