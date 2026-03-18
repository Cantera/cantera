"""Data classes common to sourcegen scaffolders."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
import re
from pathlib import Path
from sys import version_info

if version_info.minor < 11:
    from typing import Any, Iterator
    from typing_extensions import Self
else:
    from typing import Any, Iterator, Self

from ._helpers import with_unpack_iter


@dataclass
@with_unpack_iter
class Recipe:
    """
    Represents a recipe for a CLib method.

    Class holds contents of YAML header configuration.
    """

    name: str  #: name of method (without prefix)
    brief: str  #: override brief description from doxygen documentation
    what: str  #: override auto-detection of recipe type
    declaration: str  #: override auto-generated CLib declaration (custom code only)
    parameters: list  #: override auto-detected parameter doc-strings
    returns: str  #: override auto-detected return doc-string
    uses: str | list[str]  #: auxiliary C++ methods used by recipe
    wraps: str  #: signature of implemented C++ function/method
    code: str  #: custom code

    prefix: str  #: prefix used for CLib access function
    base: str  #: C++ class implementing method (if applicable)
    parents: list[str]  #: list of C++ parent classes (if applicable)
    derived: dict[str, str]  #: dictionary of C++ specialization/prefix (if applicable)

    @property
    def bases(self) -> list[str]:
        """Return all bases of a recipe."""
        return [self.base] + self.parents + list(self.derived.keys())


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
    def from_str(cls: Self, param: str, has_name: bool, doc: str = "") -> Self:
        """Generate Param from parameter string."""
        param = param.strip()
        default = None
        if "=" in param:
            param, _, default = param.partition("=")
        if has_name:  # 'param' includes both type and name
            p_type, name = param.strip().rsplit(" ", 1)
            return cls(p_type, name, doc, "", default)
        return cls(param)

    @classmethod
    def from_xml(cls: Self, param: str, has_name: bool) -> Self:
        """
        Generate Param from XML string.

        Note: Converts from doxygen style to simplified C++ whitespace notation.
        """
        for rep in [(" &", "& "), ("< ", "<"), (" >", ">"), (" *", "* ")]:
            param = param.replace(*rep)
        return cls.from_str(param.strip(), has_name=has_name)

    @classmethod
    def to_const(cls: Self, pp: Self) -> Self:
        """Create a constant version of the same parameter."""
        if pp.p_type.startswith("const "):
            return pp
        p_type = f"const {pp.p_type}"
        return cls(p_type, pp.name, pp.description, pp.direction, pp.default, pp.base)

    def short_str(self) -> str:
        """String representation of the parameter without parameter name."""
        return self.p_type

    def long_str(self) -> str:
        """String representation of the parameter with parameter name."""
        if not self.name:
            raise ValueError(f"Parameter name is undefined: {self}")
        if self.base:
            return f"{self.p_type} {self.base}::{self.name}"
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
    def from_str(cls: Self, arglist: str, has_names: bool) -> Self:
        """Generate ArgList from string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_str(arg, has_names) for arg in arglist.split(",")], spec)

    @classmethod
    def from_xml(cls: Self, arglist: str) -> Self:
        """Generate ArgList from XML string argument list."""
        arglist, spec = cls._split_arglist(arglist)
        if not arglist:
            return cls([], spec)
        return cls([Param.from_xml(arg, has_name=True) for arg in arglist.split(",")],
                   spec)

    def __len__(self) -> int:
        return len(self.params)

    def __getitem__(self, k: int) -> Param:
        return self.params[k]

    def __iter__(self) -> Iterator[Param]:
        return iter(self.params)

    def param_names_str(self) -> str:
        """String representation of the argument list parameter names."""
        args = ", ".join(par.name for par in self.params)
        return f"({args})"

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
    """Represents a function declaration.

    Func objects are a fundamental building block for code generation, where they
    hold information on C/C++ functions or their destination language equivalents. In
    addition to information used for function declarations, Func objects also contain
    annotations based on docstrings, where attribute names correspond to information
    retrieved from Doxygen tags. Future updates to ``sourcegen`` may introduce
    additional attributes, for example ``details``, ``remark`` and ``exception``, which
    can be obtained from C++ docstrings via the Doxygen XML tree.
    """

    ret_type: str  #: Return type; may include leading specifier
    name: str  #: Function name
    arglist: ArgList  #: Argument list
    brief: str = ""  #: Brief description (optional)
    wraps: Self | Param | str = None  #: Implemented C++ function/method (optional)
    returns: str = ""  #: Description of returned value (optional)
    base: str = ""  #: Qualified scope of function/method (optional)
    uses: list[Self] | None = None  #: List of auxiliary C++ methods (optional)
    deprecated: str | None = None  #: Deprecation message (if applicable)

    @classmethod
    def from_str(cls: Self, func: str, has_names: bool, brief: str = "") -> Self:
        """Generate Func from declaration string of a function declaration."""
        func = func.split("\n")[-1].rstrip(";").strip()
        # match all characters before an opening parenthesis "(" or end of line
        name = re.findall(r".*?(?=\(|$)", func)[0]
        arglist = ArgList.from_str(func.replace(name, "").strip(), has_names)
        r_type = ""
        if " " in name:
            r_type, name = name.rsplit(" ", 1)
        return cls(r_type, name, arglist, brief, None, "", "", [])

    def invocation(self) -> str:
        """Return a string representation of calling the function in an expression."""
        return (f"{self.name}{self.arglist.param_names_str()}")

    def declaration(self) -> str:
        """Return a string representation of the function without semicolon."""
        return (f"{self.ret_type} {self.name}{self.arglist.long_str()}").strip()

    @classmethod
    def from_recipe(cls: Self, recipe: Recipe) -> Self:
        """Generate annotated Func from recipe."""
        func = Func.from_str(recipe.declaration, has_names=True)

        uses = [Func.from_str(uu, has_names=False) for uu in recipe.uses]

        doc_args = []
        for arg in func.arglist:
            if arg.name not in recipe.parameters:
                msg = f"Recipe does not specify doc-string for parameter {arg.name!r}."
                raise KeyError(msg)
            doc_args.append(Param.from_str(arg.long_str(), has_name=True,
                                           doc=recipe.parameters[arg.name]))

        return cls(func.ret_type, func.name, ArgList(doc_args), recipe.brief,
                   recipe.code, recipe.returns, "", uses)

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
class HeaderFile:
    """Represents information about a parsed C header file."""

    path: Path  #: output folder
    funcs: list[Func]  #: list of functions to be scaffolded

    prefix: str = ""  #: prefix used for CLib function names
    base: str = ""  #: base class of C++ methods (if applicable)
    parents: list[str] | None = None  #: list of C++ parent class(es)
    derived: dict[str, str] | None = None  #: dictionary with alternative prefixes
    recipes: list[Recipe] | None = None  #: list of header recipes read from YAML
    docstring: list[str] | None = None  #: lines representing docstring of YAML file

    def output_name(self, suffix: str = "") -> Path:
        """
        Return output path for the generated file.

        The name of the generated file is based on the YAML configuration file
        name, where a different suffix is used. For example,
        ``<myfile>.yaml`` becomes ``<myfile>.cpp`` if the suffix is ``.cpp``.
        """
        auto, sep, suffix = suffix.partition(".")
        ret = self.path.parent / (self.path.name + auto)
        return ret.with_suffix(sep + suffix)
