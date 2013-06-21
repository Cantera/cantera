import SI, energy, mass

units = ['J__kg_K', 'kJ__kg_K', 'cal__g_K']

J__kg_K = SI.joule/(SI.kilogram * SI.kelvin)
kJ__kg_K = 1000.0*J__kg_K
cal__g_K = energy.calorie/(mass.gram * SI.kelvin)
