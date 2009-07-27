#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: energy.py,v $
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

from SI import joule

#
# Definitions of common energy units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993


Btu = 1055 * joule
erg = 1e-7 * joule
foot_pound = 1.356 * joule
horse_power_hour = 2.685e6 * joule

calorie = 4.186 * joule
Calorie = 1000 * calorie
kilowatt_hour = 3.6e6 * joule

electron_volt = 1.602e-19 * joule

#
# End of file
