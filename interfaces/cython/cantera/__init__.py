# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import os
import sys
from pathlib import Path
import warnings

from ._cantera import *  # import Cython module members
from ._cantera import (
    __git_commit__, __sundials_version__, __version__,
)

from .composite import (
    DustyGas, Interface, Quantity, Solution, SolutionArray,
)
from .data import list_data_files
from .liquidvapor import (
    CarbonDioxide, Heptane, Hfc134a, Hydrogen, Methane, Nitrogen, Oxygen, Water,
)
from .onedim import (
    BurnerFlame, CounterflowDiffusionFlame, CounterflowPremixedFlame,
    CounterflowTwinPremixedFlame, FlameBase, FreeFlame, ImpingingJet,
)
from .utils import add_module_directory, import_phases


import cantera.interrupts  # Helps with standalone packaging (PyInstaller etc.)

warnings.filterwarnings('default', module='cantera')
add_directory(Path(__file__).parent / "data")
add_directory('.')  # Move current working directory to the front of the path

# Python interpreter used for converting mechanisms
if 'PYTHON_CMD' not in os.environ:
    os.environ['PYTHON_CMD'] = sys.executable

del os
del sys
del Path
del warnings
