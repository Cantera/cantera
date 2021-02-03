from pathlib import Path
import cantera

from .test_composite import *
from .test_convert import *
from .test_equilibrium import *
from .test_func1 import *
from .test_kinetics import *
from .test_mixture import *
from .test_onedim import *
from .test_purefluid import *
from .test_reaction import *
from .test_reactor import *
from .test_thermo import *
from .test_transport import *
from .test_utils import *

cantera.add_directory(Path(__file__).parent / "data")
cantera.add_directory(Path(__file__).parents[1] / "examples" / "surface_chemistry")
