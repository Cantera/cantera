from .length import meter, centimeter, foot, inch

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
