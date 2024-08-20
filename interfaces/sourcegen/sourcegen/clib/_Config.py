"""Configuration data class used by CLib source generator."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the CLibSourceGenerator class"""

    ret_type_crosswalk = {
        "void": "int",
        "int": "int",
        "size_t": "int",
        "double": "double",
        "shared_ptr<T>": "int",
        "string": "char*",
        "vector<double>": "double[]",
        "vector<int>": "int[]",
    }

    prop_type_crosswalk = {
        "int": "int",
        "size_t": "int",
        "double": "double",
        "const string&": "const char*",
        "shared_ptr<T>": "int",
        "const shared_ptr<T>": "int",
        "const vector<double>&": 'const double*',
    }

    cabinets: List[str]

    @classmethod
    def from_parsed(cls, *, cabinets=None) -> 'Config':
        """Ensure that configurations are correct."""
        return cls(cabinets or [])
