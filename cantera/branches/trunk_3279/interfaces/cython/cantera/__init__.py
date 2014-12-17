from ._cantera import *
from ._cantera import __version__, _have_sundials
from .liquidvapor import *
from .onedim import *
from .utils import *

import os as _os
import sys as _sys
add_directory(_os.path.join(_os.path.dirname(__file__), 'data'))

# Python interpreter used for converting mechanisms
if 'PYTHON_CMD' not in _os.environ:
    _os.environ['PYTHON_CMD'] = _sys.executable
