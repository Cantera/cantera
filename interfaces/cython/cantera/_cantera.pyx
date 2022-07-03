import sys
import importlib
import importlib.abc

# Chooses the right init function
# See https://stackoverflow.com/a/52714500
class CythonPackageMetaPathFinder(importlib.abc.MetaPathFinder):
    def __init__(self, name_filter):
        super(CythonPackageMetaPathFinder, self).__init__()
        self.name_filter =  name_filter

    def find_spec(self, fullname, path, target=None):
        if fullname.startswith(self.name_filter):
            # use this extension-file but PyInit-function of another module:
            loader = importlib.machinery.ExtensionFileLoader(fullname, __file__)
            return importlib.util.spec_from_loader(fullname, loader)


# injecting custom finder/loaders into sys.meta_path:
def bootstrap_cython_submodules():
    sys.meta_path.append(CythonPackageMetaPathFinder('cantera.'))

bootstrap_cython_submodules()

from ._onedim import *
from ._utils import *
from .base import *
from .delegator import *
from .func1 import *
from .kinetics import *
from .mixture import *
from .reaction import *
from .reactionpath import *
from .reactor import *
from .speciesthermo import *
from .thermo import *
from .transport import *
from .units import *
from .yamlwriter import *
from .constants import *
