# This example solves a plug flow reactor problem, where the chemistry
# is surface chemistry. The specific problem simulated is the partial
# oxidation of methane over a platinum catalyst in a packed bed
# reactor.

from Cantera import *
from Cantera.Reactor import *
from Cantera import rxnpath
import math
import sys


#######################################################################

# unit conversion factors to SI

cm = 0.01
minute = 60.0



#######################################################################
#
# Input Parameters
#
#######################################################################

tc = 800.0                        # Temperature in Celsius

length = 0.3 * cm                 # Catalyst bed length
area = 1.0 * cm * cm              # Catalyst bed area
cat_area_per_vol = 1000.0 / cm    # Catalyst particle surface area
                                  #   per unit volume
velocity = 40.0 * cm / minute     # gas velocity
porosity = 0.3                    # Catalyst bed porosity

# input file containing the surface reaction mechanism
cti_file = 'methane_pox_on_pt.cti'

# The PFR will be simulated by a chain of 'NReactors' stirred
# reactors.
NReactors = 201
#
#  Decreased the time step by a factor of 100 to help convergence
#  12/28/2009 HKM
#
# dt = 1.0
dt = 0.01

#####################################################################


t = tc + 273.15  # convert to Kelvin

# import the gas model
gas = importPhase(cti_file,'gas')

# set the initial conditions
gas.set(T = t, P = OneAtm, X = 'CH4:1, O2:1.5, AR:0.1')
rho0 = gas.density()
nsp = gas.nSpecies()
g_names = gas.speciesNames()

# import the surface model
surf = importInterface(cti_file,'Pt_surf', [gas])
surf.setTemperature(t)
s_names = surf.speciesNames()
nsurf = surf.nSpecies()

rlen = length/(NReactors-1)
rvol = area * rlen * porosity

names = gas.speciesNames()

f = open('surf_pfr_output.csv','w')
writeCSV(f, ['Distance (mm)', 'T (C)', 'P (atm)'] + g_names  + s_names)

# catalyst area in one reactor
cat_area = cat_area_per_vol*rvol

mass_flow_rate = velocity * rho0 * area

# The plug flow reactor is represented by a linear chain of
# zero-dimensional reactors. The gas at the inlet to the first one has
# the specified inlet composition, and for all others the inlet
# composition is fixed at the composition of the reactor immediately
# upstream. Since in a PFR model there is no diffusion, the upstream
# reactors are not affected by any downstream reactors, and therefore
# the problem may be solved by simply marching from the first to last
# reactor, integrating each one to steady state.

for n in range(NReactors):

    # create a new reactor
    r = Reactor(contents = gas, energy = 'off', volume = rvol)

    # create a reservoir to represent the reactor immediately
    # upstream. Note that the gas object is set already to the
    # state of the upstream reactor
    upstream = Reservoir(gas, name = 'upstream')

    # create a reservoir for the reactor to exhaust into. The
    # composition of this reservoir is irrelevant.
    downstream = Reservoir(gas, name = 'downstream')

    # use a 'Wall' object to implement the reacting surface in the
    # reactor.  Since walls have to be installed between two
    # reactors/reserviors, we'll install it between the upstream
    # reservoir and the reactor.  The area is set to the desired
    # catalyst area in the reactor, and surface reactions are
    # included only on the side facing the reactor.
    w = Wall(left = upstream, right = r, A = cat_area, kinetics = [None, surf])
    # We need a valve between the reactor and the downstream reservoir.
    # This will determine the pressure in the reactor. Set Kv large
    # enough that the pressure difference is very small.
    v = Valve(upstream = r, downstream = downstream, Kv = 3.0e-6)

    # The mass flow rate into the reactor will be fixed by using a
    # MassFlowController object.
    m = MassFlowController(upstream = upstream,
                           downstream = r, mdot = mass_flow_rate)

    sim = ReactorNet([upstream, r, downstream])

    # set relative and absolute tolerances on the simulation
    sim.setTolerances(rtol = 1.0e-4, atol = 1.0e-11)

    time = 0
    while 1 > 0:
        time = time + dt
        sim.advance(time)

        # check whether surface coverages are in steady
        # state. This will be the case if the creation and
        # destruction rates for a surface (but not gas) species
        # are equal.
        alldone = 1

        # Note: netProduction = creation - destruction. By
        # supplying the surface object as an argument, only the
        # values for the surface species are returned by these
        # methods
        sdot = surf.netProductionRates(surf)
        cdot = surf.creationRates(surf)
        ddot = surf.destructionRates(surf)
        for ks in range(nsurf):
            ratio = sdot[ks]/(cdot[ks] + ddot[ks])
            if ratio < 0.0: ratio = -ratio
            if ratio > 1.0e-9 or time < 10*dt:
                alldone = 0

        if alldone: break

    # set the gas object state to that of this reactor, in
    # preparation for the simulation of the next reactor
    # downstream, where this object will set the inlet conditions
    gas = r.contents()

    dist = n*rlen * 1.0e3   # distance in mm

    # write the gas mole fractions and surface coverages
    # vs. distance
    writeCSV(f, [dist, r.temperature() - 273.15,
                 r.pressure()/OneAtm] + list(gas.moleFractions())
             + list(surf.coverages()))

f.close()

# make a reaction path diagram tracing carbon. This diagram will show
# the pathways by the carbon entering the bed in methane is convered
# into CO and CO2. The diagram will be specifically for the exit of
# the bed; if the pathways are desired at some interior point, then
# put this statement inside the above loop.
#
# To process this diagram, give the command on the command line
# after running this script:
# dot -Tps < carbon_pathways.dot > carbon_pathways.ps
# This will generate the diagram in Postscript.

element = 'C'
rxnpath.write(surf, element, 'carbon_pathways.dot')
