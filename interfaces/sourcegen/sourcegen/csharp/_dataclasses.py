# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass

from .._helpers import with_unpack_iter
from .._dataclasses import Func


@dataclass(frozen=True)
@with_unpack_iter
class CsFunc(Func):
    """Represents a C# interop method"""

    is_handle_release_func: bool
    handle_class_name: str | None

    def unsafe(self):
        """Identify pointers within argument lists."""
        return any(p.p_type.endswith("*") for p in self.arglist)
