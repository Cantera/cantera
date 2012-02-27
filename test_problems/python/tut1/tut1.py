#################################
print """

   Tutorial 1:              Getting started

"""
##################################


# Put this statement at the top of each Python script to import the
# most commonly-used parts of Cantera:

from Cantera import *

# The first thing you need is an object representing some phase of
# matter. We'll create here a gas mixture:
gas1 = GRI30()

# To view the state of the mixture, just print it:
print gas1

# You should see something like this:
#
#        temperature             300  K
#           pressure          101325  Pa
#            density        0.081889  kg/m^3
#   mean mol. weight         2.01588  amu

#                           1 kg            1 kmol
#                        -----------      ------------
#           enthalpy         26470.1        5.336e+04     J
#    internal energy    -1.21088e+06       -2.441e+06     J
#            entropy           64914        1.309e+05     J/K
#     Gibbs function    -1.94477e+07        -3.92e+07     J
#  heat capacity c_p         14311.8        2.885e+04     J/K
#  heat capacity c_v         10187.3        2.054e+04     J/K

#                            X                 Y
#                      -------------     ------------
#                 H2   1.000000e+00     1.000000e+00
#                  H   0.000000e+00     0.000000e+00
#                  O   0.000000e+00     0.000000e+00
#                 O2   0.000000e+00     0.000000e+00
#                 OH   0.000000e+00     0.000000e+00
#                H2O   0.000000e+00     0.000000e+00
#                HO2   0.000000e+00     0.000000e+00
#               H2O2   0.000000e+00     0.000000e+00
#                  C   0.000000e+00     0.000000e+00
#                 CH   0.000000e+00     0.000000e+00
#                CH2   0.000000e+00     0.000000e+00
#             CH2(S)   0.000000e+00     0.000000e+00
#                CH3   0.000000e+00     0.000000e+00
#                CH4   0.000000e+00     0.000000e+00
#                 CO   0.000000e+00     0.000000e+00
#                CO2   0.000000e+00     0.000000e+00
#                HCO   0.000000e+00     0.000000e+00
#               CH2O   0.000000e+00     0.000000e+00
#              CH2OH   0.000000e+00     0.000000e+00
#               CH3O   0.000000e+00     0.000000e+00
#              CH3OH   0.000000e+00     0.000000e+00
#                C2H   0.000000e+00     0.000000e+00
#               C2H2   0.000000e+00     0.000000e+00
#               C2H3   0.000000e+00     0.000000e+00
#               C2H4   0.000000e+00     0.000000e+00
#               C2H5   0.000000e+00     0.000000e+00
#               C2H6   0.000000e+00     0.000000e+00
#               HCCO   0.000000e+00     0.000000e+00
#              CH2CO   0.000000e+00     0.000000e+00
#              HCCOH   0.000000e+00     0.000000e+00
#                  N   0.000000e+00     0.000000e+00
#                 NH   0.000000e+00     0.000000e+00
#                NH2   0.000000e+00     0.000000e+00
#                NH3   0.000000e+00     0.000000e+00
#                NNH   0.000000e+00     0.000000e+00
#                 NO   0.000000e+00     0.000000e+00
#                NO2   0.000000e+00     0.000000e+00
#                N2O   0.000000e+00     0.000000e+00
#                HNO   0.000000e+00     0.000000e+00
#                 CN   0.000000e+00     0.000000e+00
#                HCN   0.000000e+00     0.000000e+00
#               H2CN   0.000000e+00     0.000000e+00
#               HCNN   0.000000e+00     0.000000e+00
#               HCNO   0.000000e+00     0.000000e+00
#               HOCN   0.000000e+00     0.000000e+00
#               HNCO   0.000000e+00     0.000000e+00
#                NCO   0.000000e+00     0.000000e+00
#                 N2   0.000000e+00     0.000000e+00
#                 AR   0.000000e+00     0.000000e+00
#               C3H7   0.000000e+00     0.000000e+00
#               C3H8   0.000000e+00     0.000000e+00
#             CH2CHO   0.000000e+00     0.000000e+00
#             CH3CHO   0.000000e+00     0.000000e+00
#
# What you have just done is to create an object ("gas1") that
# implements GRI-Mech 3.0, the 53-species, 325-reaction natural gas
# combustion mechanism developed by Gregory P. Smith, David M. Golden,
# Michael Frenklach, Nigel W. Moriarty, Boris Eiteneer, Mikhail
# Goldenberg, C. Thomas Bowman, Ronald K. Hanson, Soonho Song, William
# C. Gardiner, Jr., Vitali V. Lissianski, and Zhiwei Qin. See
# http://www.me.berkeley.edu/gri_mech/ for more information.
#
# The object created by GI30() has properties you would expect for a gas
# mixture - it has a temperature, a pressure, species mole and mass
# fractions, etc. As we'll soon see, it has many more properties.
#
# The summary of the state of 'gas1' printed above shows that new
# objects created by function GRI30() start out with a temperature of
# 300 K, a pressure of 1 atm, and have a composition that consists of
# only one species, in this case hydrogen. There is nothing special
# about H2 - it just happens to be the first species listed in the
# input file defining GRI-Mech 3.0 that the 'GRI30' function reads. In
# general, whichever species is listed first will initially have a
# mole fraction of 1.0, and all of the others will be zero.


#  Setting the state
#  -----------------

# The state of the object can easily be changed. For example,

gas1.setTemperature(1200)

# sets the temperature to 1200 K. (Cantera always uses SI units.)
# After this statement,

print gas1

# results in:
#
#       temperature                 1200  K
#       pressure                  405300  Pa
#       density                 0.081896  kg/m^3
#       mean mol. weight         2.01594  amu
#
#                           X                 Y
#                     -------------     ------------
#                H2   1.000000e+000     1.000000e+000
#                (other species not shown)
#
# Notice that the temperature has been changed as requested, but the
# pressure has changed too. The density and composition have
# not.
#
# When setting properties individually, some convention needs to be
# adopted to specify which other properties are held constant. This is
# because thermodynamics requires that *two* properties (not one) in
# addition to composition information be specified to fix the
# intensive state of a substance (or mixture).
#
# Cantera adopts the following convention: only one of the set
# (temperature, density, mass fractions) is altered by setting any
# single property. This means that:
#
# a) Setting the temperature is done holding density and
#    composition  fixed. (The pressure changes.)

# b) Setting the pressure is done holding temperature and
#    composition fixed. (The density changes.)
#
# c) Setting the composition is done holding temperature
#    and density fixed. (The pressure changes).
#

# Instead of using a method like 'setTemperature' to set one property,
# you can use a single method 'set' to set any property or combination
# of properties:

gas1.set(Temperature = 900.0, Pressure = 1.e5)

# This statement sets both temperature and pressure at the same
# time. Any number of property/value pairs can be specified in a
# call to 'set'. For example, the following sets the mole fractions
# too:

gas1.set(Temperature = 900.0, Pressure = 1.e5,
         MoleFractions = 'CH4:1,O2:2,N2:7.52')

# The 'set' function also accepts abbreviated property names:

gas1.set(T = 900.0, P = 1.0e5, X = 'CH4:1,O2:2,N2:7.52')

# Either version results in:
print gas1

#        temperature             900  K
#           pressure          100000  Pa
#            density        0.369279  kg/m^3
#   mean mol. weight         27.6332  amu

#                           1 kg            1 kmol
#                        -----------      ------------
#           enthalpy          455660        1.259e+07     J
#    internal energy          184862        5.108e+06     J
#            entropy         8529.31        2.357e+05     J/K
#     Gibbs function    -7.22072e+06       -1.995e+08     J
#  heat capacity c_p          1304.4        3.604e+04     J/K
#  heat capacity c_v         1003.52        2.773e+04     J/K

#                            X                 Y
#                      -------------     ------------
#                 H2   0.000000e+00     0.000000e+00
#                  H   0.000000e+00     0.000000e+00
#                  O   0.000000e+00     0.000000e+00
#                 O2   1.901141e-01     2.201487e-01
#                 OH   0.000000e+00     0.000000e+00
#                H2O   0.000000e+00     0.000000e+00
#                HO2   0.000000e+00     0.000000e+00
#               H2O2   0.000000e+00     0.000000e+00
#                  C   0.000000e+00     0.000000e+00
#                 CH   0.000000e+00     0.000000e+00
#                CH2   0.000000e+00     0.000000e+00
#             CH2(S)   0.000000e+00     0.000000e+00
#                CH3   0.000000e+00     0.000000e+00
#                CH4   9.505703e-02     5.518632e-02
#                 CO   0.000000e+00     0.000000e+00
#                CO2   0.000000e+00     0.000000e+00
#                HCO   0.000000e+00     0.000000e+00
#               CH2O   0.000000e+00     0.000000e+00
#              CH2OH   0.000000e+00     0.000000e+00
#               CH3O   0.000000e+00     0.000000e+00
#              CH3OH   0.000000e+00     0.000000e+00
#                C2H   0.000000e+00     0.000000e+00
#               C2H2   0.000000e+00     0.000000e+00
#               C2H3   0.000000e+00     0.000000e+00
#               C2H4   0.000000e+00     0.000000e+00
#               C2H5   0.000000e+00     0.000000e+00
#               C2H6   0.000000e+00     0.000000e+00
#               HCCO   0.000000e+00     0.000000e+00
#              CH2CO   0.000000e+00     0.000000e+00
#              HCCOH   0.000000e+00     0.000000e+00
#                  N   0.000000e+00     0.000000e+00
#                 NH   0.000000e+00     0.000000e+00
#                NH2   0.000000e+00     0.000000e+00
#                NH3   0.000000e+00     0.000000e+00
#                NNH   0.000000e+00     0.000000e+00
#                 NO   0.000000e+00     0.000000e+00
#                NO2   0.000000e+00     0.000000e+00
#                N2O   0.000000e+00     0.000000e+00
#                HNO   0.000000e+00     0.000000e+00
#                 CN   0.000000e+00     0.000000e+00
#                HCN   0.000000e+00     0.000000e+00
#               H2CN   0.000000e+00     0.000000e+00
#               HCNN   0.000000e+00     0.000000e+00
#               HCNO   0.000000e+00     0.000000e+00
#               HOCN   0.000000e+00     0.000000e+00
#               HNCO   0.000000e+00     0.000000e+00
#                NCO   0.000000e+00     0.000000e+00
#                 N2   7.148289e-01     7.246650e-01
#                 AR   0.000000e+00     0.000000e+00
#               C3H7   0.000000e+00     0.000000e+00
#               C3H8   0.000000e+00     0.000000e+00
#             CH2CHO   0.000000e+00     0.000000e+00
#             CH3CHO   0.000000e+00     0.000000e+00


# Other properties may also be set using 'set', including some that
# can only be set in combination with others.  The following property
# pairs may be set: (Enthalpy, Pressure), (IntEnergy, Volume),
# (Entropy, Volume), (Entropy, Pressure). In each case, the values of
# the extensive properties must be entered *per unit mass*.

# Setting the enthalpy and pressure:
gas1.set(Enthalpy = 2*gas1.enthalpy_mass(), Pressure = 2*OneAtm)

# This sets gas1 to a state with P = 2 atm, and a specific enthalpy
# twice its previous value.

# Note that the abbreviations T, P, H, U, S, V can also be used with
# the 'set' method.

# The composition above was specified using a string. The format is a
# comma-separated list of <species name>:<relative mole numbers>
# pairs. The mole numbers will be normalized to produce the mole
# fractions, and therefore they are 'relative' mole numbers.  Mass
# fractions can be set in this way too by changing 'X' to 'Y' in the
# above statement.

# The composition can also be set using an array, which must have the
# same size as the number of species. For example, to set all 53 mole
# fractions to the same value, do this:

x = ones(53,'d');   # NumPy array of 53 ones
gas1.set(X = x)
print gas1

# To set the mass fractions to equal values:
gas1.set(Y = x)
print gas1
