# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import sourcegen

sourcegen.generate_source(*sys.argv[1:])
