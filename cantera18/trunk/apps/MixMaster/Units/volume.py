#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: volume.py,v $
# Revision 1.1.1.1  2003/04/14 17:57:49  dggoodwin
# Initial import.
#
# Revision 1.1  2002/12/20 13:23:41  dgg
# first commit
#
# Revision 1.1.1.1  2000/01/21 22:59:51  dgg
# dgg Cantera
#
# Revision 1.1  1999/11/27 17:27:35  dgg
# initial import to Cantera repository
#
# Revision 1.1  1999/11/25 19:50:58  aivazis
# Original source
#

from length import meter, centimeter, foot, inch

#
# Definitions of common volume units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993

cubic_meter = meter**3
cubic_centimeter = centimeter**3
cubic_foot = foot**3
cubic_inch = inch**3

liter = 1000 * cubic_centimeter

us_fluid_ounce = 231./128 * cubic_inch
us_pint = 16 * us_fluid_ounce
us_fluid_quart = 2 * us_pint
us_fluid_gallon = 4 * us_fluid_quart

#
# End of file
