# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from pathlib import Path
import re

from typing import List, Any, Tuple
try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml

from ._helpers import with_unpack_iter


BlockMap = yaml.comments.CommentedMap

@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """Represents a function parameter"""

    p_type: str
    name: str = ""
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
            return Param(*parts, default)
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
    def _split_arglist(arglist: str) -> Tuple[str, str]:
        arglist = arglist.strip()
        suffix = arglist[arglist.rfind(")") + 1:]
        arglist = re.findall(re.compile(r'(?<=\().*(?=\))'), arglist)[0]
        return arglist, suffix.strip()

    @staticmethod
    def from_str(arglist: str) -> 'ArgList':
        """Generate ArgList from string argument list"""
        arglist, suffix = ArgList._split_arglist(arglist)
        if not arglist:
            return ArgList([], suffix)
        return ArgList([Param.from_str(_) for _ in arglist.split(",")], suffix)

    @staticmethod
    def from_xml(arglist: str) -> 'ArgList':
        """Generate ArgList from XML string argument list"""
        arglist, suffix = ArgList._split_arglist(arglist)
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
    """Represents a function parsed from a C header file"""

    annotations: str
    ret_type: str
    name: str
    params: List[Param]


@dataclass(frozen=True)
@with_unpack_iter
class AnnotatedFunc(Func):
    """Represents a function annotated with doxygen info."""

    implements: str
    relates: str
    cxx_type: str
    cxx_name: str
    cxx_arglist: str
    cxx_anchorfile: str
    cxx_anchor: str

    @classmethod
    def to_yaml(cls, representer, func):
        out = BlockMap([
            ('annotations', yaml.scalarstring.PreservedScalarString(func.annotations)),
            ('ret_type', func.ret_type),
            ('name', func.name),
            ('params', f"({', '.join([_.p_type for _ in func.params])})"),
            ('implements', func.implements),
            ('relates', func.relates),
            ('cxx_type', Param.from_xml(func.cxx_type).short_str()),
            ('cxx_name', func.cxx_name),
            ('cxx_arglist', ArgList.from_xml(func.cxx_arglist).long_str()),
            ('cxx_anchorfile', func.cxx_anchorfile),
            ('cxx_anchor', func.cxx_anchor),
        ])
        return representer.represent_dict(out)


@dataclass(frozen=True)
@with_unpack_iter
class HeaderFile:
    """Represents information about a parsed C header file"""

    path: Path
    funcs: List[Func]
