# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
import re
from pathlib import Path
from typing import List, Any, Tuple, Iterator

from ._helpers import with_unpack_iter


@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """Represents a function parameter"""

    p_type: str
    name: str = ""

    @classmethod
    def from_str(cls, param: str) -> 'Param':
        """Generate Param from string parameter"""
        parts = param.strip().rsplit(" ", 1)
        if len(parts) == 2 and parts[0] not in ["const", "virtual", "static"]:
            return cls(*parts)
        return cls(param)

    def long_str(self) -> str:
        """String representation of the parameter with parameter name."""
        if not self.name:
            raise ValueError(f"Parameter name is undefined: {self}")
        return f"{self.p_type} {self.name}"


@dataclass(frozen=True)
class ArgList:
    """Represents a function argument list"""

    params: List[Param]
    spec: str = ""  #: trailing Specification (example: `const`)

    @staticmethod
    def _split_arglist(arglist: str) -> Tuple[str, str]:
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

    def __getitem__(self, k):
        return self.params[k]

    def __iter__(self) -> "Iterator[Param]":
        return iter(self.params)

    def long_str(self) -> str:
        """String representation of the argument list with parameter names."""
        args = ', '.join([par.long_str() for par in self.params])
        return f"({args}) {self.spec}".strip()


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """Represents a function parsed from a C header file."""

    ret_type: str  # may include leading specifier
    name: str
    arglist: ArgList

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
class HeaderFile:
    """Represents information about a parsed C header file"""

    path: Path
    funcs: List[Func]
