import gc
import numpy as np
import pytest
from pytest import approx

import cantera as ct
phase = ct.Solution('oxygen-plasma.yaml', 'isotropic-electron-energy-plasma',transport_model=None)
phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
print(phase.elastic_power_loss)
