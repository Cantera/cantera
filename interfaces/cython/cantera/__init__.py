# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

from ._cantera import *
from ._cantera import __version__, __sundials_version__, __git_commit__
from .composite import *
from .liquidvapor import *
from .onedim import *
from .utils import *

import os
import sys
add_directory(os.path.join(os.path.dirname(__file__), 'data'))

# Python interpreter used for converting mechanisms
if 'PYTHON_CMD' not in os.environ:
    os.environ['PYTHON_CMD'] = sys.executable
