# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
import re
from pathlib import Path
from typing import List, Any

from ._helpers import with_unpack_iter


@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """Represents a function parameter"""

    p_type: str
    name: str = ""
    description: str = ""
    direction: str = ""
    default: Any = None

    @staticmethod
    def from_str(param: str) -> 'Param':
        """Generate Param from string parameter"""
        default = None
        if "=" in param:
            default = param[param.rfind("=") + 1:]
            param = param[:param.rfind("=")]
        parts = param.strip().rsplit(" ", 1)
        if len(parts) == 2 and parts[0] not in ["const", "virtual", "static"]:
            return Param(*parts, "", "", default)
        return Param(param)

    @staticmethod
    def from_xml(param: str) -> 'Param':
        replacements = [(" &amp;", "& "), ("&lt; ", "<"), (" &gt;", ">")]
        for rep in replacements:
            param = param.replace(*rep)
        ret = Param.from_str(param.strip())
        return ret

    def short_str(self) -> str:
        """Return a short string representation of the parameter"""
        return f"{self.p_type}"

    def long_str(self) -> str:
        """Return a short string representation of the parameter"""
        if not self.name:
            raise ValueError(f"Parameter name is undefined: {self}")
        if self.default is not None:
            return f"{self.p_type} {self.name}={self.default}"
        return f"{self.p_type} {self.name}"


@dataclass(frozen=True)
@with_unpack_iter
class ArgList:
    """Represents a function argument list"""

    params: List[Param]
    suffix: str = ""

    @staticmethod
    def from_str(arglist: str) -> 'ArgList':
        """Generate ArgList from string argument list"""
        arglist, suffix = split_arglist(arglist)
        if not arglist:
            return ArgList([], suffix)
        return ArgList([Param.from_str(_) for _ in arglist.split(",")], suffix)

    @staticmethod
    def from_xml(arglist: str) -> 'ArgList':
        """Generate ArgList from XML string argument list"""
        arglist, suffix = split_arglist(arglist)
        if not arglist:
            return ArgList([], suffix)
        return ArgList([Param.from_xml(_) for _ in arglist.split(",")], suffix)

    def short_str(self) -> str:
        """Return a short string representation of the argument list"""
        return f"({', '.join([_.short_str() for _ in self.params])})"

    def long_str(self) -> str:
        """Return a short string representation of the argument list"""
        return f"({', '.join([_.long_str() for _ in self.params])}) {self.suffix}".strip()

    def n_optional(self):
        """Return the number of optional arguments"""
        return sum([1 for _ in self.params if _.default is not None])


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """Represents a function parsed from a C header file."""

    annotations: str
    ret_type: str
    name: str
    params: List[Param]


@dataclass(frozen=True)
@with_unpack_iter
class Recipe:
    """Represents a recipe for a CLib method."""

    prefix: str
    name: str
    base: str
    uses: List[str]
    implements: str
    relates: str = ""
    what: str = ""


@dataclass(frozen=True)
@with_unpack_iter
class HeaderFile:
    """Represents information about a parsed C header file"""

    path: Path
    funcs: List[Func]
    recipes: List[Recipe] = None


def split_arglist(arglist: str) -> tuple:
    """Split C++ argument list into text within parentheses and suffix."""
    arglist = arglist.strip()
    suffix = arglist[arglist.rfind(")") + 1:]
    arglist = re.findall(re.compile(r'(?<=\().*(?=\))'), arglist)[0]
    return arglist, suffix.strip()
