# An equilibrium example with charged species in the gas phase
# and multiple condensed phases.

# Note: This example runs fine on Mac and linux platforms, but
# encounters some convergence difficulties under Windows. The reasons
# for this are not currently known.

from Cantera import *

# create objects representing the gas phase and the condensed
# phases. The gas is a mixture of multiple species, and the condensed
# phases are all modeled as incompressible stoichiometric
# substances. See file KOH.cti for more information.
phases = importPhases('KOH.cti',
                      ['K_solid',
                       'K_liquid', 'KOH_a', 'KOH_b',
                       'KOH_liquid', 'K2O2_solid',
                       'K2O_solid', 'KO2_solid',
                       'ice', 'liquid_water','KOH_plasma'])

# create the Mixture object from the list of phases
mix = Mixture(phases)

# open the output file and write the column headings
f = open('equil_koh.csv','w')
writeCSV(f,['T']+mix.speciesNames())

# loop over temperature
for n in range(100):
    t = 350.0 + 50.0*n
    print 'T = ',t
    mix.set(T= t, P = OneAtm, Moles="K:1.03, H2:2.12, O2:0.9")

    # set the mixture to a state of chemical equilibrium holding
    # temperature and pressure fixed
    # mix.equilibrate("TP",maxsteps=10000,loglevel=1)
    mix.vcs_equilibrate("TP",printLvl=0,maxsteps=10000,loglevel=0)

    # write out the moles of each species
    writeCSV(f,[t]+ list(mix.speciesMoles()))

# close the output file
f.close()
