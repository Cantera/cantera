"""
This example solves a plug flow reactor problem, where the chemistry is
surface chemistry. The specific problem simulated is the partial oxidation of
methane over a platinum catalyst in a packed bed reactor.

Requires: cantera >= 2.5.0
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

output_filename = 'surf_pfr_output.csv'

# The PFR will be simulated by a chain of 'NReactors' stirred reactors.
NReactors = 201
dt = 1.0

#####################################################################

t = tc + 273.15  # convert to Kelvin

# import the gas model and set the initial conditions
gas = ct.Solution(yaml_file, 'gas')
gas.TPX = t, ct.one_atm, 'CH4:1, O2:1.5, AR:0.1'

# import the surface model
surf = ct.Interface(yaml_file, 'Pt_surf', [gas])
surf.TP = t, ct.one_atm

rlen = length/(NReactors-1)
rvol = area * rlen * porosity

# catalyst area in one reactor
cat_area = cat_area_per_vol * rvol

mass_flow_rate = velocity * gas.density * area

# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.

TDY = gas.TDY
cov = surf.coverages

print('    distance       X_CH4        X_H2        X_CO')

# create a new reactor
gas.TDY = TDY
r = ct.IdealGasReactor(gas, energy='off')
r.volume = rvol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas, name='downstream')

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.
rsurf = ct.ReactorSurface(surf, r, A=cat_area)

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(r, downstream, master=m, K=1e-5)

sim = ct.ReactorNet([r])
sim.max_err_test_fails = 12

# set relative and absolute tolerances on the simulation
sim.rtol = 1.0e-9
sim.atol = 1.0e-21

output_data = []

for n in range(NReactors):
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r.thermo.TDY
    upstream.syncState()
    sim.reinitialize()
    sim.advance_to_steady_state()
    dist = n * rlen * 1.0e3  # distance in mm

    if n % 10 == 0:
        print('  {0:10f}  {1:10f}  {2:10f}  {3:10f}'.format(
            dist, *gas['CH4', 'H2', 'CO'].X))

    # write the gas mole fractions and surface coverages vs. distance
    output_data.append(
        [dist, r.T - 273.15, r.thermo.P/ct.one_atm] + list(gas.X)
        + list(surf.coverages)
    )

with open(output_filename, 'w', newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Distance (mm)', 'T (C)', 'P (atm)'] +
                    gas.species_names + surf.species_names)
    writer.writerows(output_data)

print("Results saved to '{0}'".format(output_filename))
