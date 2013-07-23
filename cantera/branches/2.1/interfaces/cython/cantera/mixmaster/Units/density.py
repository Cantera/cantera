import math
from . import SI

#
# Definitions of common density units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993
units = ['kg__m3', 'g__m3', 'g__cm3']

kg__m3 = SI.kilogram/(SI.meter * SI.meter * SI.meter)
g__m3 = 1.e-3*kg__m3
g__cm3 = 1.e3*kg__m3
