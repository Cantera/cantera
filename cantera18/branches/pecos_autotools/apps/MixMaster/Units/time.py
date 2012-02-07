#!/bin/env python
#
#--------------------------------------------------------------------------
#
# $License$
#
#--------------------------------------------------------------------------

# $Log: time.py,v $
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

from SI import second

#
# Definitions of common time units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993

minute = 60 * second
hour = 60 * minute
day = 24 * hour
year = 365.25 * day

#
# End of file
