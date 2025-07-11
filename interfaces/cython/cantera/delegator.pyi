# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from typing import Callable

from cantera.reaction import ExtensibleRateData

def extension(
    *, name: str, data: ExtensibleRateData | None = None
) -> Callable[[type], type]: ...
