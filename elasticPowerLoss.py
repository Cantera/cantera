import gc
import numpy as np
import pytest
from pytest import approx

import cantera as ct
plasma = ct.Solution('example_data/oxygen-plasma-itikawa.yaml',
                     'isotropic-electron-energy-plasma',
                      transport_model=None)
print("plasma temp: ",plasma.Te)
P = ct.one_atm * 0.01
T = 400.0
plasma.TPX = T, P, 'O2:1.0, e:5e-3, O2+:5e-3'
print("plasma temp: ",plasma.Te)
plasma.mean_electron_energy = 10 # [eV]
print("plasma temp: ",plasma.Te)
phase = ct.Solution('oxygen-plasma.yaml', 'isotropic-electron-energy-plasma',transport_model=None)
phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
print(phase.elastic_power_loss)

'''
plasma temp:  38681.7270718336
plasma temp:  38681.7270718336
plasma temp:  77363.4541436672
10269608971.625109

plasma temp:  0.0
plasma temp:  0.0
plasma temp:  1479.566489519614
'''
