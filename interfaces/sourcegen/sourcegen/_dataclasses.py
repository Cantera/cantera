"""Data classes common to sourcegen scaffolders."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
import re
from pathlib import Path
from typing import Any, Iterator
from typing_extensions import Self

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
    base: str = ""  #: Base (optional). Only used if param represents a member variable

    @classmethod
    def from_str(cls: Self, param: str, doc: str = "") -> Self:
        """Generate Param from parameter string."""
        param = param.strip()
        default = None
        if "=" in param:
            param, _, default = param.partition("=")
        parts = param.strip().rsplit(" ", 1)
        if len(parts) == 2 and parts[0] not in ["const", "virtual", "static"]:
            if "@param" not in doc:
                return cls(*parts, "", "", default)
            items = doc.split()
            if items[1] != parts[1]:
                msg = f"Documented variable {items[1]!r} does not match {parts[1]!r}"
                raise ValueError(msg)
            direction = items[0].split("[")[1].split("]")[0] if "[" in items[0] else ""
            return cls(*parts, " ".join(items[2:]), direction, default)
        return cls(param)

    @classmethod
    def from_xml(cls: Self, param: str) -> Self:
        """
        Generate Param from XML string.

        Note: Converts from doxygen style to simplified C++ whitespace notation.
        """
        for rep in [(" &", "& "), ("< ", "<"), (" >", ">"), (" *", "* ")]:
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
        regex = re.compile(r"(?<=\().*(?=\))", flags=re.DOTALL)
        arglist = re.findall(regex, arglist)[0]
        return arglist, spec

    @classmethod
    def from_str(cls: Self, arglist: str) -> Self:
        """Generate ArgList from string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_str(arg) for arg in arglist.split(",")], spec)

    @classmethod
    def from_xml(cls: Self, arglist: str) -> Self:
        """Generate ArgList from XML string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_xml(arg) for arg in arglist.split(",")], spec)

    def __len__(self) -> int:
        return len(self.params)

    def __getitem__(self, k: int) -> Param:
        return self.params[k]

    def __iter__(self) -> Iterator[Param]:
        return iter(self.params)

    def short_str(self) -> str:
        """String representation of the argument list without parameter names."""
        args = ", ".join(par.short_str() for par in self.params)
        return f"({args}) {self.spec}".strip()

    def long_str(self) -> str:
        """String representation of the argument list with parameter names."""
        args = ", ".join([par.long_str() for par in self.params])
        return f"({args}) {self.spec}".strip()


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """Represents a function declaration in a C/C++ header file."""

    ret_type: str  #: Return type; may include leading specifier
    name: str  #: Function name
    arglist: ArgList  #: Argument list

    @classmethod
    def from_str(cls: Self, func: str) -> Self:
        """Generate Func from declaration string of a function."""
        func = func.rstrip(";").strip()
        # match all characters before an opening parenthesis "(" or end of line
        name = re.findall(r".*?(?=\(|$)", func)[0]
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

    brief: str = ""  #: Brief description (optional)
    implements: Self = None  #: Implemented C++ function/method (optional)
    returns: str = ""  #: Description of returned value (optional)
    base: str = ""  #: Qualified scope of function/method (optional)
    uses: list[Self] | None = None  #: List of auxiliary C++ methods (optional)

    @classmethod
    def from_str(cls: Self, func: str, brief: str = "") -> Self:
        """Generate annotated CFunc from header block of a function."""
        lines = func.split("\n")
        func = Func.from_str(lines[-1])
        if len(lines) == 1:
            return cls(*func, brief, None, "", "", [])
        returns = ""
        doc_args = {p.name: p for p in func.arglist}
        for ix, line in enumerate(lines[:-1]):
            line = line.strip().lstrip("*").strip()
            if ix == 1 and not brief:
                brief = line
            elif line.startswith("@param"):
                # match parameter name
                keys = [k for k in doc_args.keys() if line.split()[1] == k]
                if len(keys) == 1:
                    key = keys[0]
                    doc_args[key] = Param.from_str(doc_args[key].long_str(), line)
            elif line.startswith("@returns"):
                returns = line.lstrip("@returns").strip()
        args = ArgList(list(doc_args.values()))
        return cls(func.ret_type, func.name, args, brief, None, returns, "", [])

    def short_declaration(self) -> str:
        """Return a short string representation."""
        if self.arglist is None:
            ret = (f"{self.name}").strip()
        else:
            ret = (f"{self.name}{self.arglist.short_str()}").strip()
        if self.base:
            return f"{self.ret_type} {self.base}::{ret}"
        return f"{self.ret_type} {ret}"

    @property
    def ret_param(self) -> Param:
        """Assemble return parameter."""
        return Param(self.ret_type, "", self.returns)


@dataclass
@with_unpack_iter
class Recipe:
    """
    Represents a recipe for a CLib method.

    Class holds contents of YAML header configuration.
    """

    name: str  #: name of method (without prefix)
    implements: str  #: signature of implemented C++ function/method
    uses: str | list[str]  #: auxiliary C++ methods used by recipe
    what: str  #: override auto-detection of recipe type
    brief: str  #: override brief description from doxygen documentation
    code: str  #: custom code to override autogenerated code (stub: to be implemented)

    prefix: str  #: prefix used for CLib access function
    base: str  #: C++ class implementing method (if applicable)
    parents: list[str]  #: list of C++ parent classes (if applicable)
    derived: list[str]  #: list of C++ specializations (if applicable)


@dataclass
class HeaderFile:
    """Represents information about a parsed C header file."""

    path: Path  #: output folder
    funcs: list[Func]  #: list of functions to be scaffolded

    prefix: str = ""  #: prefix used for CLib function names
    base: str = ""  #: base class of C++ methods (if applicable)
    parents: list[str] | None = None  #: list of C++ parent class(es)
    derived: list[str] | None = None  #: list of C++ specialization(s)
    recipes: list[Recipe] | None = None  #: list of header recipes read from YAML
    docstring: list[str] | None = None  #: lines representing docstring of YAML file

    def output_name(self, suffix: str = "") -> Path:
        """
        Return output path for the generated file.

        The name of the auto-generated file is based on the YAML configuration file
        name, where ``_auto`` is stripped and a different suffix is used. For example,
        ``<myfile>_auto.yaml`` becomes ``<myfile>3.cpp`` if the suffix is ``3.cpp``.
        """
        auto, sep, suffix = suffix.partition(".")
        ret = self.path.parent / self.path.name.replace("_auto", auto)
        return ret.with_suffix(sep + suffix)
