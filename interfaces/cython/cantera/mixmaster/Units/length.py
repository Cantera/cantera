from .SI import meter
from .SI import nano, milli, centi, kilo

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
