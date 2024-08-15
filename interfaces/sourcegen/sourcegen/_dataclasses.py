# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from pathlib import Path
from typing import List
import re
from ._helpers import with_unpack_iter


@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """Represents a function parameter"""

    p_type: str
    name: str


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """Represents a function parsed from a C header file"""

    ret_type: str
    name: str
    params: List[Param]
    comments: str


@dataclass(frozen=True)
@with_unpack_iter
class HeaderFile:
    """Represents information about a parsed C header file"""

    path: Path
    funcs: List[Func]
