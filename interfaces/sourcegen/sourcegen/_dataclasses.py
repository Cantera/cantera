# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List
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
    name: str


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
            ('cxx_type', func.cxx_type),
            ('cxx_name', func.cxx_name),
            ('cxx_arglist', func.cxx_arglist),
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
