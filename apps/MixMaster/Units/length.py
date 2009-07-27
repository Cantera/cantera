#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: length.py,v $
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

from SI import meter
from SI import nano, milli, centi, kilo

#
# Definitions of common length units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993


nanometer = nano * meter
millimeter = milli * meter
centimeter = centi * meter
kilometer = kilo * meter

inch = 2.540 * centimeter
foot = 12 * inch
yard = 3 * foot
mile = 5280 * foot

fathom = 6 * foot
nautical_mile = 6076 * foot

angstrom = 1e-10 * meter
fermi = 1e-15 * meter
light_year = 9.460e12 * kilometer
parsec = 3.084e13 * kilometer

#
# End of file
