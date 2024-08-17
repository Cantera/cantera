# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the YamlSourceGenerator class"""

    bases: List[str]

    @classmethod
    def from_parsed(cls, *, bases: List[str]=None) -> 'Config':
        """Ensure that configurations are correct."""
        return cls(bases or [])
