# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Cython modules built as their own extensions (e.g. ``jacobians``, enhancement #241)
# reference shared symbols such as ``pyCanteraError`` that are defined in the merged
# ``_cantera`` extension. Load ``_cantera`` (and the submodule extensions it pulls in)
# with RTLD_GLOBAL so those symbols resolve across extension boundaries, then restore the
# previous flags. POSIX-only; on Windows symbols are exported explicitly instead.
import sys as _sys
import os as _os
_old_dlopenflags = _sys.getdlopenflags() if hasattr(_sys, "getdlopenflags") else None
if _old_dlopenflags is not None:
    _sys.setdlopenflags(_old_dlopenflags | _os.RTLD_GLOBAL)
try:
    from ._cantera import *
finally:
    if _old_dlopenflags is not None:
        _sys.setdlopenflags(_old_dlopenflags)
from ._utils import __version__, __sundials_version__, __git_commit__
from .composite import *
from .liquidvapor import *
from .onedim import *
from .data import *
import cantera.interrupts  # Helps with standalone packaging (PyInstaller etc.)

from pathlib import Path
import warnings

warnings.filterwarnings('default', module='cantera')
add_data_directory(Path(__file__).parent / "data")
add_data_directory('.')  # Move current working directory to the front of the path

del warnings, Path, np, os
