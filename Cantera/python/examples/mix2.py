# Mixing two streams with reaction. This is the same as mix1.py,
# except that a source of H atoms is added to ignite the fuel/air
# mixture. Once ignited, the flow of H atoms is stopped.

import math
from Cantera import *
from Cantera.Reactor import Reactor, Reservoir, MassFlowController, Valve


# Use air for stream a. Note that the Air() function does not set the
# composition correctly; thus, we need to explicitly set the
# composition to that of air.
gas_a = Air()
gas_a.setState_TPX(300.0, OneAtm, 'O2:0.21, N2:0.78, AR:0.01')
rho_a = gas_a.density()


# Use GRI-Mech 3.0 for stream b (methane) and for the mixer. If it is
# desired to have a pure mixer, with no chemistry, use instead a
# reaction mechanism for gas_b that has no reactions.
gas_b = GRI30()
gas_b.setState_TPX(300.0, OneAtm, 'CH4:1')
rho_b = gas_b.density()


# Create reservoirs for the two inlet streams and for the outlet
# stream.  The upsteam reservoirs could be replaced by reactors, which
# might themselves be connected to reactors further upstream. The
# outlet reservoir could be replaced with a reactor with no outlet, if
# it is desired to integrate the composition leaving the mixer in
# time, or by an arbitrary network of downstream reactors.
res_a = Reservoir(gas_a)
res_b = Reservoir(gas_b)
downstream   = Reservoir(gas_b)


# Create a reactor for the mixer. A reactor is required instead of a
# reservoir, since the state will change with time if the inlet mass
# flow rates change or if there is chemistry occurring.
mixer = Reactor(gas_b)


# create two mass flow controllers connecting the upstream reservoirs
# to the mixer, and set their mass flow rates to values corresponding
# to stoichiometric combustion.
mfc1 = MassFlowController(res_a, mixer)
mfc1.setMassFlowRate(rho_a*2.5/0.21)

mfc2 = MassFlowController(res_b, mixer)
mfc2.setMassFlowRate(rho_b*1.0)


# connect the mixer to the downstream reservoir with a valve.
outlet = Valve(mixer, downstream)
outlet.setValveCoeff(1.0)


# add an igniter to ignite the mixture. The 'igniter' consists of a
# stream of pure H.
gas_c = IdealGasMix('h2o2.xml')
gas_c.setState_TPX(300.0, OneAtm, 'H:1')
igniter = Reactor(gas_c)

mfc3 = MassFlowController(igniter, mixer)
mfc3.setMassFlowRate(0.05)



# Since the mixer is a reactor, we need to integrate in time to reach
# steady state. A few residence times should be enough.
t = 0.0
for n in range(30):
    tres = mixer.mass()/(mfc1.massFlowRate() + mfc2.massFlowRate())
    tnow = t
    t += 0.5*tres
    mixer.advance(t)
    
    # if ignited, turn the igniter off.
    # We also need to restart the integration in this case.
    if mixer.temperature() > 1200.0:
        mfc3.setMassFlowRate(0.0)
        mixer.setInitialTime(t)
        
    print '%14.5g %14.5g %14.5g  %14.5g  %14.5g' % (t, mixer.temperature(),
                                                    mixer.enthalpy_mass(),
                                                    mixer.pressure(),
                                                    mixer.massFraction('CH4'))
gas_b.setState_TPY(mixer.temperature(), mixer.pressure(), mixer.massFractions())

# view the state of the gas in the mixer
gas_b.setState_TPY(mixer.temperature(), mixer.pressure(),
                   mixer.massFractions())
print gas_b    
