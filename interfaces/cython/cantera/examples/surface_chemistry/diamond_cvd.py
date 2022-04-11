"""
A CVD example simulating growth of a diamond film

This example computes the growth rate of a diamond film according to a
simplified version of a particular published growth mechanism (see file
diamond.yaml for details). Only the surface coverage equations are solved here;
the gas composition is fixed. (For an example of coupled gas-phase and
surface, see catalytic_combustion.py.) Atomic hydrogen plays an important
role in diamond CVD, and this example computes the growth rate and surface
coverages as a function of [H] at the surface for fixed temperature and [CH3].

Requires: cantera >= 2.6.0, pandas >= 0.25.0, matplotlib >= 2.0
Keywords: surface chemistry, kinetics
"""

import csv
import cantera as ct

print('\n******  CVD Diamond Example  ******\n')

# import the model for the diamond (100) surface and the adjacent bulk phases
d = ct.Interface("diamond.yaml", "diamond_100")
g = d.adjacent["gas"]
dbulk = d.adjacent["diamond"]

mw = dbulk.molecular_weights[0]

t = 1200.0
x = g.X
p = 20.0 * ct.one_atm / 760.0  # 20 Torr
g.TP = t, p

ih = g.species_index('H')

xh0 = x[ih]

with open('diamond.csv', 'w', newline='') as f:
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

    print('H concentration, growth rate, and surface coverages '
          'written to file diamond.csv')

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    data = pd.read_csv('diamond.csv')

    data.plot(x="H mole Fraction", y="Growth Rate (microns/hour)", legend=False)
    plt.xlabel('H Mole Fraction')
    plt.ylabel('Growth Rate (microns/hr)')
    plt.show()

    names = [name for name in data.columns if not name.startswith(('H mole', 'Growth'))]
    data.plot(x='H mole Fraction', y=names, legend=True)
    plt.xlabel('H Mole Fraction')
    plt.ylabel('Coverage')
    plt.show()
except ImportError:
    print("Install matplotlib and pandas to plot the outputs")
