from ._cantera import *
from ._cantera import __version__, __sundials_version__
from .liquidvapor import *
from .onedim import *
from .utils import *

import os
import sys
add_directory(os.path.join(os.path.dirname(__file__), 'data'))

# Python interpreter used for converting mechanisms
if 'PYTHON_CMD' not in os.environ:
    os.environ['PYTHON_CMD'] = sys.executable
