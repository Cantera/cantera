# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from typing import Union

from .._helpers import with_unpack_iter
from .._dataclasses import Func


@dataclass(frozen=True)
@with_unpack_iter
class DoxyFunc(Func):
    """Represents a Doxygen function info."""

    implements: str
