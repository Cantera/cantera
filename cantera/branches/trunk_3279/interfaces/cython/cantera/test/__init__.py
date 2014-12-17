import os
import cantera

cantera.add_directory(os.path.join(os.path.dirname(__file__), 'data'))

from .test_thermo import *
from .test_purefluid import *
from .test_equilibrium import *
from .test_kinetics import *
from .test_transport import *
from .test_mixture import *
from .test_func1 import *
from .test_reactor import *
from .test_onedim import *
from .test_convert import *