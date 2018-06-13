"""
A CVD example simulating growth of a diamond film

This example computes the growth rate of a diamond film according to a
simplified version of a particular published growth mechanism (see file
diamond.cti for details). Only the surface coverage equations are solved here;
the gas composition is fixed. (For an example of coupled gas- phase and
surface, see catalytic_combustion.py.)  Atomic hydrogen plays an important
role in diamond CVD, and this example computes the growth rate and surface
coverages as a function of [H] at the surface for fixed temperature and [CH3].
"""

import csv
import cantera as ct

print('\n******  CVD Diamond Example  ******\n')

# import the models for the gas and bulk diamond
g, dbulk = ct.import_phases('diamond.cti', ['gas', 'diamond'])

# import the model for the diamond (100) surface
d = ct.Interface('diamond.cti', 'diamond_100', [g, dbulk])

ns = d.n_species
mw = dbulk.molecular_weights[0]

t = 1200.0
x = g.X
p = 20.0 * ct.one_atm / 760.0  # 20 Torr
g.TP = t, p

ih = g.species_index('H')

xh0 = x[ih]
f = open('diamond.csv', 'w')
writer = csv.writer(f)
writer.writerow(['H mole Fraction', 'Growth Rate (microns/hour)'] +
                d.species_names)

iC = d.kinetics_species_index(dbulk.species_index('C(d)'), 1)

for n in range(20):
    x[ih] /= 1.4
    g.TPX = t, p, x
    d.advance_coverages(10.0)  # integrate the coverages to steady state
    carbon_dot = d.net_production_rates[iC]
    mdot = mw * carbon_dot
    rate = mdot / dbulk.density
    writer.writerow([x[ih], rate * 1.0e6 * 3600.0] + list(d.coverages))

f.close()

print('H concentration, growth rate, and surface coverages '
      'written to file diamond.csv')
