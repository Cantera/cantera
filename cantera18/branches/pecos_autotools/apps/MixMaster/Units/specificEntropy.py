#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: specificEntropy.py,v $
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

units = ['J__kg_K', 'kJ__kg_K', 'cal__g_K']

J__kg_K = SI.joule/(SI.kilogram * SI.kelvin)
kJ__kg_K = 1000.0*J__kg_K
cal__g_K = energy.calorie/(mass.gram * SI.kelvin)

#
# End of file
