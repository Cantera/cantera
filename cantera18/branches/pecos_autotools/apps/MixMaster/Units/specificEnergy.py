#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: specificEnergy.py,v $
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

import SI, energy, mass

#
# Definitions of common energy units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993
units = ['J__kg', 'kJ__kg', 'Btu__lbm', 'cal__g', 'kcal__g', 'kcal__kg']

J__kg = SI.joule/SI.kilogram
kJ__kg = 1000.0*J__kg
Btu__lbm = energy.Btu/mass.pound
cal__g = energy.calorie/mass.gram
kcal__g = 1000.0*cal__g
kcal__kg = cal__g
#
# End of file
