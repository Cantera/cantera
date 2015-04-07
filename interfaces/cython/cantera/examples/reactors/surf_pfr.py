"""
This example solves a plug flow reactor problem, where the chemistry is
surface chemistry. The specific problem simulated is the partial oxidation of
methane over a platinum catalyst in a packed bed reactor.
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
cti_file = 'methane_pox_on_pt.cti'

output_filename = 'surf_pfr_output.csv'

# The PFR will be simulated by a chain of 'NReactors' stirred reactors.
NReactors = 201
dt = 1.0

#####################################################################

t = tc + 273.15  # convert to Kelvin

# import the gas model and set the initial conditions
gas = ct.Solution(cti_file, 'gas')
gas.TPX = t, ct.one_atm, 'CH4:1, O2:1.5, AR:0.1'

# import the surface model
surf = ct.Interface(cti_file,'Pt_surf', [gas])
surf.TP = t, ct.one_atm

rlen = length/(NReactors-1)
rvol = area * rlen * porosity

outfile = open(output_filename,'w')
writer = csv.writer(outfile)
writer.writerow(['Distance (mm)', 'T (C)', 'P (atm)'] +
                gas.species_names + surf.species_names)

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

for n in range(NReactors):
    surf.TP = TDY[0], ct.one_atm
    surf.coverages = cov

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

    # use a 'Wall' object to implement the reacting surface in the reactor.
    # Since walls have to be installed between two reactors/reserviors, we'll
    # install it between the upstream reservoir and the reactor.  The area is
    # set to the desired catalyst area in the reactor, and surface reactions
    # are included only on the side facing the reactor.
    w = ct.Wall(upstream, r, A=cat_area, kinetics=[None, surf])

    # We need a valve between the reactor and the downstream reservoir. This
    # will determine the pressure in the reactor. Set K large enough that the
    # pressure difference is very small.
    v = ct.Valve(r, downstream, K=1e-4)

    # The mass flow rate into the reactor will be fixed by using a
    # MassFlowController object.
    m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)

    sim = ct.ReactorNet([r])
    sim.max_err_test_fails = 12

    # set relative and absolute tolerances on the simulation
    sim.rtol = 1.0e-9
    sim.atol = 1.0e-21

    T_start, rho_start, Y_start = r.thermo.TDY
    cov_start = surf.coverages
    V_start = r.volume

    Tu_start, rhou_start, Yu_start = upstream.thermo.TDY

    time = 0
    all_done = False
    while not all_done:
        time += dt
        sim.advance(time)

        # check whether surface coverages are in steady state. This will be
        # the case if the creation and destruction rates for a surface (but
        # not gas) species are equal.
        all_done = True

        # Note: netProduction = creation - destruction. By supplying the
        # surface object as an argument, only the values for the surface
        # species are returned by these methods
        sdot = surf.get_net_production_rates(surf)
        cdot = surf.get_creation_rates(surf)
        ddot = surf.get_destruction_rates(surf)

        for ks in range(surf.n_species):
            ratio = abs(sdot[ks]/(cdot[ks] + ddot[ks]))
            if ratio > 1.0e-9 or time < 10*dt:
                all_done = False

    # Save the reactor and surface states, in preparation for the simulation
    # of the next reactor downstream, where this object will set the inlet
    # conditions and the initial surface coverages
    TDY = r.thermo.TDY
    cov = surf.coverages

    dist = n * rlen * 1.0e3   # distance in mm

    if not n % 10:
        print('  {0:10f}  {1:10f}  {2:10f}  {3:10f}'.format(dist, *gas['CH4','H2','CO'].X))

    # write the gas mole fractions and surface coverages vs. distance
    writer.writerow([dist, r.T - 273.15, r.thermo.P/ct.one_atm] +
                    list(gas.X) + list(surf.coverages))

outfile.close()
print("Results saved to '{0}'".format(output_filename))
