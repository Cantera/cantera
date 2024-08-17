# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from typing import Dict

from .._helpers import get_preamble, normalize_indent


@dataclass(frozen=True)
class Config:
    """Provides configuration info for the CSharpSourceGenerator class"""


    func_prolog = normalize_indent("""
        [DllImport(LibFile)]
        public static extern
    """)

    ret_type_crosswalk = {
        "const char*": "string",
        "const double*": "double[]",
        "size_t": "nuint",
        "char*": "byte*"
    }

    prop_type_crosswalk = {
        "byte*": "string",
        "double*": "double[]"
    }

    preamble = "/*\n" + get_preamble() + "*/"

    # These we load from the parsed YAML config file
    class_crosswalk: Dict[str, str]

    derived_handles: Dict[str, str]

    wrapper_classes: Dict[str, Dict[str, str]]

    @classmethod
    def from_parsed(cls, *,
                    class_crosswalk=None, derived_handles=None, wrapper_classes=None):
        return cls(class_crosswalk or {}, derived_handles or {}, wrapper_classes or {})
