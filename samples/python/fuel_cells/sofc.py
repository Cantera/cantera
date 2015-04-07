# SOFC
#
# This script implements a simple model of a solid oxide fuel
# cell. Unlike most SOFC models, however, it does not use
# semi-empirical Butler-Volmer kinetics for the charge transfer
# reactions, but uses elementary, reversible reactions obeying
# mass-action kinetics for all reactions, including charge
# transfer. As this script will demonstrate, this approach allows
# computing the OCV (it does not need to be separately specified), as
# well as polarization curves.
#
# NOTE: The parameters here, and in the input file sofc.cti, are not
# to be relied upon for a real SOFC simulation! They are meant to
# illustrate only how to do such a calculation in Cantera. While some
# of the parameters may be close to real values, others are simply set
# arbitrarily to give reasonable-looking results.

# It is recommended that you read input file sofc.cti before reading
# or running this script!

#---------------------------------------------------------------------


from Cantera import *
import math


#--------------------------------------------------
#
#   parameters
#
#--------------------------------------------------

tc = 800.0                  # T in C
temp = tc + 273.15
pres = OneAtm

# gas compositions. Change as desired.
anode_gas_X   = 'H2:0.97, H2O:0.03'
cathode_gas_X = 'O2:1.0, H2O:0.001'

# time to integrate coverage eqs. to steady state in
# 'advanceCoverages'. This should be more than enough time.
tss = 50.0

# electrolyte conductivity
sigma =  2.0        # Siemens / m

# electrolyte thickness
ethick = 5.0e-5     # m

# TPB length per unit area
TPB_length_per_area = 1.0e7   # per meter



#----------------------------------------------------
#
#   utility functions
#
#----------------------------------------------------

def showCoverages(s):
    """Print the coverages for surface s."""
    print '\n '+s.name() + '\n'
    cov = s.coverages()
    names = s.speciesNames()
    nsp = len(names)
    for n in range(nsp):
        print '%16s  %13.4g ' % (names[n], cov[n])

def equil_OCV(gas1, gas2):
    return -GasConstant*gas1.temperature()*math.log(gas1.moleFraction('O2')/
                                              gas2.moleFraction('O2'))/(4.0*Faraday)

def NewtonSolver(f, xstart, C = 0.0):
    """Solve f(x) = C by Newton iteration.
    - xstart    starting point for Newton iteration
    - C         constant
    """
    f0 = f(xstart) - C
    x0 = xstart
    dx = 1.0e-6
    xlast = 999.0
    n = 0
    while n < 200:
        ff = f(x0 + dx) - C
        dfdx = (ff - f0)/dx
        step = - f0/dfdx

        # avoid taking steps too large
        if abs(step) > 0.1:
            step = 0.1*step/abs(step)

        x0 += step
        emax = 0.00001  # 0.01 mV tolerance
        if abs(f0) < emax and n > 8:
            return x0
        xlast = x0
        f0 = f(x0) - C
        n += 1
    raise 'no root!'



#####################################################################
#
# Anode-side phases
#
#####################################################################

# import the anode-side bulk phases
gas_a, anode_bulk, oxide_a = importPhases('sofc.cti',
                                          ['gas', 'metal', 'oxide_bulk',])

# import the surfaces on the anode side
anode_surf = importInterface('sofc.cti','metal_surface',[gas_a])
oxide_surf_a = importInterface('sofc.cti','oxide_surface',[gas_a, oxide_a])

# import the anode-side triple phase boundary
tpb_a = importEdge('sofc.cti', 'tpb', [anode_bulk, anode_surf, oxide_surf_a])

anode_surf.setName('anode surface')
oxide_surf_a.setName('anode-side oxide surface')

# this function is defined to use with NewtonSolver to invert the
# current-voltage function. NewtonSolver requires a function of one
# variable, so the other objects are accessed through the global
# namespace.
def anode_curr(E):
    """Current from the anode as a function of anode
    potential relative to electrolyte"""

    # the anode-side electrolyte potential is kept at zero.
    # Therefore, the anode potential is just equal to E.
    anode_bulk.setElectricPotential(E)

    # get the species net production rates due to the anode-side TPB
    # reaction mechanism. The production rate array has the values for
    # the neighbor species in the order listed in the .cti file,
    # followed by the tpb phase. Since the first neighbor phase is the
    # bulk metal, species 0 is the electron.
    w = tpb_a.netProductionRates()

    # the sign convention is that the current is positive when
    # electrons are being delivered to the anode - i.e. it is positive
    # for fuel cell operation.
    return Faraday * w[0] * TPB_length_per_area




#####################################################################
#
# Cathode-side phases
#
#####################################################################

# Here for simplicity we are using the same phase and interface models
# for the cathode as we used for the anode. In a more realistic
# simulation, separate models would be used for the cathode, with a
# different reaction mechanism.

# import the cathode-side bulk phases
gas_c, cathode_bulk, oxide_c = importPhases('sofc.cti',
                                            ['gas',
                                             'metal',
                                             'oxide_bulk',])

# import the surfaces on the cathode side
cathode_surf = importInterface('sofc.cti','metal_surface',[gas_c])
oxide_surf_c = importInterface('sofc.cti','oxide_surface',[gas_c, oxide_c])

# import the cathode-side triple phase boundary
tpb_c = importEdge('sofc.cti', 'tpb', [cathode_bulk,
                                       cathode_surf, oxide_surf_c])

cathode_surf.setName('cathode surface')
oxide_surf_c.setName('cathode-side oxide surface')

def cathode_curr(E):
    """Current to the cathode as a function of cathode
    potential relative to electrolyte"""

    # due to ohmic losses, the cathode-side electrolyte potential is
    # non-zero. Therefore, we need to add this potential to E to get
    # the cathode potential.
    ee = E + oxide_c.electricPotential()
    cathode_bulk.setElectricPotential(ee)

    # get the species net production rates due to the cathode-side TPB
    # reaction mechanism. The production rate array has the values for
    # the neighbor species in the order listed in the .cti file,
    # followed by the tpb phase. Since the first neighbor phase is the
    # bulk metal, species 0 is the electron.
    w = tpb_c.netProductionRates()

    # the sign convention is that the current is positive when electrons
    # are being drawn from the cathode (i.e, negative production rate).
    return -Faraday * w[0] * TPB_length_per_area



# initialization

# set the gas compositions, and temperatures of all phases

gas_a.set(T = temp, P = pres, X = anode_gas_X)
gas_a.equilibrate('TP')  # needed to use equil_OCV

gas_c.set(T = temp, P = pres, X = cathode_gas_X)
gas_c.equilibrate('TP')  # needed to use equil_OCV

phases = [anode_bulk, anode_surf, oxide_surf_a, oxide_a, cathode_bulk,
          cathode_surf, oxide_surf_c, oxide_c, tpb_a, tpb_c]
for p in phases:
    p.setTemperature(temp)




# now bring the surface coverages into steady state with these gas
# compositions. Note that the coverages are held fixed at these values
# - we do NOT consider the change in coverages due to TPB
# reactions. For that, a more complex model is required. But as long
# as the thermal chemistry is fast relative to charge transfer, this
# should be an OK approximation.

for s in [anode_surf, oxide_surf_a, cathode_surf, oxide_surf_c]:
    s.advanceCoverages(tss)
    showCoverages(s)



# find open circuit potentials by solving for the E values that give
# zero current.

Ea0 = NewtonSolver(anode_curr, xstart = -0.51)
Ec0 = NewtonSolver(cathode_curr, xstart = 0.51)

print '\nocv from zero current is: ',Ec0 - Ea0
print 'OCV from thermo equil is: ',equil_OCV(gas_a, gas_c)

print 'Ea0 = ', Ea0
print 'Ec0 = ', Ec0
print

# do polarization curve for anode overpotentials from -250 mV
# (cathodic) to +250 mV (anodic)
Ea_min = Ea0 - 0.25
Ea_max = Ea0 + 0.25

file = open('sofc.csv','w')

writeCSV(file,['i (mA/cm2)','eta_a','eta_c','eta_ohmic', 'Eload'])

# vary the anode overpotential, from cathodic to anodic polarization
for n in range(100):
    Ea = Ea_min + 0.005*n

    # set the electrode potential. Note that the anode-side electrolyte
    # is held fixed at 0 V.
    anode_bulk.setElectricPotential(Ea)

    # compute the anode current
    curr = anode_curr(Ea)

    # set potential of the oxide on the cathode side to reflect
    # the ohmic drop through the electrolyte

    delta_V = curr * ethick / sigma

    # if the current is positive, negatively-charged ions are flowing
    # from the cathode to the anode. Therefore, the cathode side must be
    # more negative than the anode side.
    phi_oxide_c = -delta_V

    # note that both the bulk and the surface potentials must be set
    oxide_c.setElectricPotential(phi_oxide_c)
    oxide_surf_c.setElectricPotential(phi_oxide_c)

    # Find the value of the cathode potential relative to the
    # cathode-side electrolyte that yields the same current density
    # as the anode current density
    Ec = NewtonSolver(cathode_curr, xstart = Ec0 + 0.1, C = curr)

    cathode_bulk.setElectricPotential(phi_oxide_c + Ec);

    # write the current density, anode and cathode overpotentials,
    # ohmic overpotential, and load potential
    writeCSV(file,[0.1*curr, Ea - Ea0, Ec - Ec0, delta_V,
                   cathode_bulk.electricPotential()
                   - anode_bulk.electricPotential()])

print 'polarization curve data written to file sofc.csv'

file.close()
