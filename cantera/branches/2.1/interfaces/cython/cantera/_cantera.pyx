#cython: embedsignature=True
#distutils: language = c++

import numpy as np
cimport numpy as np
import math

from cython.operator cimport dereference as deref

from _cantera cimport *

include "utils.pyx"
include "constants.pyx"
include "func1.pyx"

include "base.pyx"
include "thermo.pyx"
include "kinetics.pyx"
include "transport.pyx"
include "composite.pyx"

include "mixture.pyx"
include "reactor.pyx"
include "onedim.pyx"
include "reactionpath.pyx"
