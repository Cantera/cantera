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
