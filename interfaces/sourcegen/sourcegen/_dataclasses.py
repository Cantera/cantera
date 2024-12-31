"""Data classes common to sourcegen scaffolders."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
import re
from pathlib import Path
from typing import Any, Iterator

from ._helpers import with_unpack_iter


@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """Class representing a function parameter."""

    p_type: str  #: Parameter type
    name: str = ""  #: Parameter name; may be empty if used for return argument
    description: str = ""  #: Parameter description (optional annotation)
    direction: str = ""  #: Direction of parameter (optional annotation)
    default: Any = None  #: Default value (optional)

    @classmethod
    def from_str(cls, param: str) -> 'Param':
        """Generate Param from parameter string."""
        param = param.strip()
        default = None
        if "=" in param:
            default = param[param.rfind("=")+1:]
            param = param[:param.rfind("=")]
        parts = param.strip().rsplit(" ", 1)
        if len(parts) == 2 and parts[0] not in ["const", "virtual", "static"]:
            return cls(*parts, "", "", default)
        return cls(param)

    @classmethod
    def from_xml(cls, param: str) -> 'Param':
        """Generate Param from XML string."""
        param = param.strip()
        replacements = [(" &amp;", "& "), ("&lt; ", "<"), (" &gt;", ">")]
        for rep in replacements:
            param = param.replace(*rep)
        return cls.from_str(param.strip())

    def short_str(self) -> str:
        """String representation of the parameter without parameter name."""
        return self.p_type

    def long_str(self) -> str:
        """String representation of the parameter with parameter name."""
        if not self.name:
            raise ValueError(f"Parameter name is undefined: {self}")
        return f"{self.p_type} {self.name}"


@dataclass(frozen=True)
class ArgList:
    """Represents a function argument list."""

    params: list[Param]  #: List of function parameters
    spec: str = ""  #: Trailing specification (example: `const`)

    @staticmethod
    def _split_arglist(arglist: str) -> tuple[str, str]:
        """Split string into text within parentheses and trailing specification."""
        arglist = arglist.strip()
        if not arglist:
            return "", ""
        spec = arglist[arglist.rfind(")") + 1:]
        # match text within parentheses
        regex = re.compile(r'(?<=\().*(?=\))', flags=re.DOTALL)
        arglist = re.findall(regex, arglist)[0]
        return arglist, spec

    @classmethod
    def from_str(cls, arglist: str) -> 'ArgList':
        """Generate ArgList from string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_str(arg) for arg in arglist.split(",")], spec)

    @classmethod
    def from_xml(cls, arglist: str) -> 'ArgList':
        """Generate ArgList from XML string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_xml(arg) for arg in arglist.split(",")], spec)

    def __len__(self):
        return len(self.params)

    def __getitem__(self, k):
        return self.params[k]

    def __iter__(self) -> "Iterator[Param]":
        return iter(self.params)

    def short_str(self) -> str:
        """String representation of the argument list without parameter names."""
        args = ', '.join([par.short_str() for par in self.params])
        return f"({args}) {self.spec}".strip()

    def long_str(self) -> str:
        """String representation of the argument list with parameter names."""
        args = ', '.join([par.long_str() for par in self.params])
        return f"({args}) {self.spec}".strip()


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """Represents a function declaration in a C/C++ header file."""

    ret_type: str  #: Return type; may include leading specifier
    name: str  #: Function name
    arglist: ArgList  #: Argument list

    @classmethod
    def from_str(cls, func: str) -> 'Func':
        """Generate Func from declaration string of a function."""
        func = func.strip()
        # match all characters before an opening parenthesis '(' or end of line
        name = re.findall(r'.*?(?=\(|$)', func)[0]
        arglist = ArgList.from_str(func.replace(name, "").strip())
        r_type = ""
        if " " in name:
            r_type, name = name.rsplit(" ", 1)
        return cls(r_type, name, arglist)

    def declaration(self) -> str:
        """Return a string representation of the function without semicolon."""
        return (f"{self.ret_type} {self.name}{self.arglist.long_str()}").strip()


@dataclass(frozen=True)
@with_unpack_iter
class CFunc(Func):
    """Represents an annotated function declaration in a C/C++ header file."""

    brief: str=""  #: Brief description (optional)
    implements: 'Func'=None  #: Implemented C++ function/method (optional)
    returns: str=""  #: Description of returned value (optional)
    base: str=""  #: Qualified scope of function/method (optional)
    relates: list[str]=None  #: List of related C++ methods (optional)

    def short_declaration(self) -> str:
        """Return a short string representation."""
        ret = (f"{self.name}{self.arglist.short_str()}").strip()
        if self.base:
            return f"{self.ret_type} {self.base}::{ret}"
        return f"{self.ret_type} {ret}"

    @property
    def ret_param(self):
        """Assemble return parameter."""
        return Param(self.ret_type, "", self.returns)


@dataclass
@with_unpack_iter
class Recipe:
    """
    Represents a recipe for a CLib method.

    Class holds contents of YAML header configuration.
    """

    prefix: str  #: Prefix used for CLib access function
    name: str  #: Name of method (without prefix)
    base: str  #: C++ class implementing method
    parents: list[str]  #: List of C++ parent class(es)
    derived: list[str]  #: List of C++ specializations
    uses: list[str]  #: List of referenced CLib cabinets
    implements: str  #: Signature of implemented method
    relates: list[str]=None  #: Methods used to retrieve instances of managed objects
    size_fcn: str=""  #: Method used to check array size
    what: str = ""  #: Non-empty for special methods: "constructor", "destructor"


@dataclass
# @with_unpack_iter
class HeaderFile:
    """Represents information about a parsed C header file"""

    path: Path
    funcs: list[Func]
    recipes: list[Recipe] = None

    def output_name(self, auto="3", suffix=""):
        """Return updated path."""
        ret = self.path.parent / self.path.name.replace("_auto", auto)
        return ret.with_suffix(suffix)
