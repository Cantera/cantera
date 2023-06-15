# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# This code to set the application registry has to come before any code that wants to
# use the registry is imported. In particular, it has to come before any of our code
# that uses units!
from pint import UnitRegistry, set_application_registry
cantera_units_registry = UnitRegistry()
set_application_registry(cantera_units_registry)

# Now we can import our code
from .solution import *
