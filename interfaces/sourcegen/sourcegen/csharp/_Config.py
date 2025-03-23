# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from sys import version_info

if version_info.minor < 11:
    from typing_extensions import Self
else:
    from typing import Self


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the CSharpSourceGenerator class."""


    ret_type_crosswalk = {
        "const char*": "string",
        "const double*": "double[]",
        "const int*": "int[]",
        "size_t": "nuint",
        "char*": "byte*"
    }

    prop_type_crosswalk = {
        "byte*": "string",
        "double*": "double[]"
    }

    # These we load from the parsed YAML config file
    class_crosswalk: dict[str, str]

    class_accessors: dict[str, str]

    derived_handles: dict[str, str]

    wrapper_classes: dict[str, dict[str, str]]

    @classmethod
    def from_parsed(cls: Self, *,
                    class_crosswalk: dict[str, str] | None = None,
                    class_accessors: dict[str, str] | None = None,
                    derived_handles: dict[str, str] | None = None,
                    wrapper_classes: dict[str, dict[str, str]] | None = None):
        return cls(class_crosswalk or {}, class_accessors or {},
                   derived_handles or {}, wrapper_classes or {})
