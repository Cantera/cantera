"""
This example solves a plug flow reactor problem with surface chemistry. The specific
problem simulated is the partial oxidation of methane over a platinum catalyst in a
packed bed reactor. This example solves the DAE system directly, using the FlowReactor
class and the SUNDIALS IDA solver, in contrast to the approximation as a chain of
steady-state WSRs used in surf_pfr_chain.py.

Requires: cantera >= 3.0.0
Keywords: catalysis, reactor network, surface chemistry, plug flow reactor,
          packed bed reactor
"""

import csv

import cantera as ct

# unit conversion factors to SI
cm = 0.01
minute = 60.0

#######################################################################
# Input Parameters
#######################################################################

tc = 800.0  # Temperature in Celsius
length = 0.3 * cm  # Catalyst bed length
area = 1.0 * cm**2  # Catalyst bed area
cat_area_per_vol = 1000.0 / cm  # Catalyst particle surface area per unit volume
velocity = 40.0 * cm / minute  # gas velocity
porosity = 0.3  # Catalyst bed porosity

# input file containing the surface reaction mechanism
yaml_file = 'methane_pox_on_pt.yaml'

output_filename = 'surf_pfr2_output.csv'

#####################################################################

t = tc + 273.15  # convert to Kelvin

# import the model and set the initial conditions
surf = ct.Interface(yaml_file, 'Pt_surf')
surf.TP = t, ct.one_atm
gas = surf.adjacent['gas']
gas.TPX = t, ct.one_atm, 'CH4:1, O2:1.5, AR:0.1'

mass_flow_rate = velocity * gas.density * area * porosity

# create a new reactor
r = ct.FlowReactor(gas)
r.area = area
r.surface_area_to_volume_ratio = cat_area_per_vol * porosity
r.mass_flow_rate = mass_flow_rate
r.energy_enabled = False

# Add the reacting surface to the reactor
rsurf = ct.ReactorSurface(surf, r)

sim = ct.ReactorNet([r])

output_data = []
n = 0
print('    distance       X_CH4        X_H2        X_CO')
print('  {:10f}  {:10f}  {:10f}  {:10f}'.format(
      0, *r.thermo['CH4', 'H2', 'CO'].X))

while sim.distance < length:
    dist = sim.distance * 1e3  # convert to mm
    sim.step()

    if n % 100 == 0 or (dist > 1 and n % 10 == 0):
        print('  {:10f}  {:10f}  {:10f}  {:10f}'.format(
              dist, *r.thermo['CH4', 'H2', 'CO'].X))
    n += 1

    # write the gas mole fractions and surface coverages vs. distance
    output_data.append(
        [dist, r.T - 273.15, r.thermo.P / ct.one_atm]
        + list(r.thermo.X)  # use r.thermo.X not gas.X
        + list(rsurf.kinetics.coverages)  # use rsurf.kinetics.coverages not surf.coverages
    )

with open(output_filename, 'w', newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Distance (mm)', 'T (C)', 'P (atm)'] +
                    gas.species_names + surf.species_names)
    writer.writerows(output_data)

print("Results saved to '{0}'".format(output_filename))
