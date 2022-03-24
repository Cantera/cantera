# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from ._cantera import *
from ._cantera import __version__, __sundials_version__, __git_commit__
from .composite import *
from .liquidvapor import *
from .onedim import *
from .utils import *
import cantera.interrupts  # Helps with standalone packaging (PyInstaller etc.)

import os
import sys
from pathlib import Path
import warnings

warnings.filterwarnings('default', module='cantera')
add_directory(Path(__file__).parent / "data")
add_directory('.')  # Move current working directory to the front of the path

# Python interpreter used for converting mechanisms
if 'PYTHON_CMD' not in os.environ:
    os.environ['PYTHON_CMD'] = sys.executable
