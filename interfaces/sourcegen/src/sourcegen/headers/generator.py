"""Resolve recipes based on specification file contents."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import logging
from dataclasses import dataclass

from ruamel import yaml
from jinja2 import Environment, BaseLoader

from .tagfiles import TagFileParser
from ..dataclasses import HeaderFile, Param, ArgList, Func, Recipe
from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the HeaderGenerator class"""

    ret_type_crosswalk: dict[str, str]  #: Return type crosswalks

    par_type_crosswalk: dict[str, str]  #: Parameter type crosswalks


class HeaderGenerator:
    """Convert recipes to header file information."""

    _clib_bases: list[str] = None  #: list of bases provided via YAML configurations

    def __init__(self, config: dict, templates: dict, bases: list[str]) -> None:
        self._config = Config(**config)
        self._templates = templates
        self._doxygen_tags = None
        self._clib_bases = bases

    def resolve_tags(self, headers_files: list[HeaderFile], root: str) -> None:
        """Resolve recipe information based on Doxygen tags."""
        def get_bases() -> tuple[list[str], list[str]]:
            bases = set()
            classes = set()
            for headers in headers_files:
                bases.add(headers.base)
                for recipe in headers.recipes:
                    classes.update(recipe.bases)
            return sorted(bases), sorted(classes)

        self._clib_bases, classes = get_bases()
        self._doxygen_tags = TagFileParser(root, classes)

        for headers in headers_files:
            msg = f"  resolving recipes in {headers.path.name!r}:"
            _LOGGER.info(msg)
            c_funcs = []
            for recipe in headers.recipes:
                c_funcs.append(self.resolve_recipe(recipe))
            headers.funcs = c_funcs

    def resolve_recipe(self, recipe: Recipe) -> Func:
        """Build CLib header from recipe and Doxygen annotations."""
        def merge_params(
                wraps: str, cxx_member: Func | Param
            ) -> tuple[list[Param], Func]:
            """Create preliminary CLib argument list."""
            obj_handle = []
            if "::" in wraps:
                # If class method, add handle as first parameter
                what = wraps.split("::")[0]
                obj_handle.append(
                    Param("int", "handle", f"Handle to queried {what} object."))
            if isinstance(cxx_member, Param):
                if recipe.what.endswith("setter"):
                    return obj_handle + [cxx_member], cxx_member
                return obj_handle, cxx_member

            if "(" not in wraps:
                return obj_handle + cxx_member.arglist.params, cxx_member

            # Signature may skip C++ default parameters
            args_short = Func.from_str(wraps).arglist
            if len(args_short) < len(cxx_member.arglist):
                cxx_arglist = ArgList(cxx_member.arglist[:len(args_short)])
                cxx_member = Func(cxx_member.ret_type, cxx_member.name,
                                  cxx_arglist, cxx_member.brief, cxx_member.wraps,
                                  cxx_member.returns, cxx_member.base, cxx_member.uses)

            return obj_handle + cxx_member.arglist.params, cxx_member

        func_name = f"{recipe.prefix}_{recipe.name}"
        reserved = ["cabinetSize", "resetStorage"]
        if recipe.name in reserved:
            loader = Environment(loader=BaseLoader)
            msg = f"   generating {func_name!r} -> {recipe.what}"
            _LOGGER.debug(msg)
            header = loader.from_string(
                self._templates[f"clib-reserved-{recipe.name}"]
                ).render(base=recipe.base, prefix=recipe.prefix)
            reader = yaml.YAML(typ="safe")
            header = reader.load(header)
            for key, value in header.items():
                recipe.__setattr__(key, value)
            return Func.from_recipe(recipe)

        if recipe.code:
            # Custom code
            return Func.from_recipe(recipe)

        # Ensure that all functions/methods referenced in recipe are detected correctly
        bases = recipe.bases
        if not recipe.wraps:
            recipe.wraps = self._doxygen_tags.detect(recipe.name, bases)
        elif recipe.base and "::" not in recipe.wraps:
            parts = list(recipe.wraps.partition("("))
            parts[0] = self._doxygen_tags.detect(parts[0], bases)
            recipe.wraps = "".join(parts)
        recipe.uses = [self._doxygen_tags.detect(uu.split("(")[0], bases, False)
                       for uu in recipe.uses]

        cxx_member = None
        ret_param = Param("void")
        args = []
        brief = ""

        if recipe.wraps:
            cxx_member = self._doxygen_tags.cxx_member(
                recipe.wraps, recipe.what.endswith("setter"))

            if cxx_member.base in recipe.derived:
                # Use alternative prefix for class specialization
                recipe.prefix = recipe.derived.get(cxx_member.base, recipe.prefix)
                func_name = f"{recipe.prefix}_{recipe.name}"

            msg = f"   generating {func_name!r} -> {recipe.wraps}"
            _LOGGER.debug(msg)

            if isinstance(cxx_member, Func):
                # Convert C++ return type to format suitable for crosswalk:
                # Incompatible return parameters are buffered and appended to back
                ret_param, buffer_params = self._ret_crosswalk(
                    cxx_member.ret_type, recipe.derived)
                par_list, cxx_member = merge_params(recipe.wraps, cxx_member)
                prop_params = self._prop_crosswalk(par_list)
                args = prop_params + buffer_params
                brief = cxx_member.brief
            elif recipe.what == "variable-setter":
                ret_param = Param("int32_t")
                par_list, cxx_member = merge_params(recipe.wraps, cxx_member)
                args = self._prop_crosswalk(par_list)
                if args[-1].p_type == "char*":
                    args[-1] = Param.to_const(args[-1])
                brief = cxx_member.description
            else:
                # Variable getter
                prop_params, cxx_member = merge_params(recipe.wraps, cxx_member)
                ret_param, buffer_params = self._ret_crosswalk(
                    cxx_member.p_type, recipe.derived)
                args = prop_params + buffer_params
                brief = cxx_member.description

        if recipe.what and cxx_member:
            # Recipe type and corresponding C++ function are known
            pass

        elif cxx_member and isinstance(cxx_member, Param):
            # Recipe represents a variable getter/setter
            recipe.what = "variable-getter"

        elif cxx_member:
            # Autodetection of CLib function purpose ("what")
            cxx_arglen = len(cxx_member.arglist)
            if not cxx_member.base:
                if (cxx_member.name.startswith("new") and
                    any(base in cxx_member.ret_type
                        for base in [recipe.base] + list(recipe.derived.keys()))):
                    recipe.what = "constructor"
                else:
                    recipe.what = "function"
            elif "void" not in cxx_member.ret_type and cxx_arglen == 0:
                recipe.what = "getter"
            elif "void" in cxx_member.ret_type and cxx_arglen == 1:
                p_type = cxx_member.arglist[0].p_type
                if cxx_member.name.startswith("get"):
                    recipe.what = "getter"
                elif "*" in p_type and not p_type.startswith("const"):
                    recipe.what = "getter"  # getter assigns to existing array
                else:
                    recipe.what = "setter"
            elif any(recipe.wraps.startswith(base) for base in recipe.bases):
                recipe.what = "method"
            else:
                msg = f"Unable to auto-detect function type for recipe {recipe.name!r}."
                _LOGGER.critical(msg)
                exit(1)

        elif recipe.what == "noop":
            # No operation
            msg = f"   generating {func_name!r} -> no-operation"
            _LOGGER.debug(msg)
            args = [Param("int32_t", "handle", f"Handle to {recipe.base} object.")]
            brief = "No operation."
            ret_param = Param("int32_t", "", "Always zero.")

        elif recipe.name == "new":
            # Default constructor
            recipe.what = "constructor"
            msg = f"   generating {func_name!r} -> default constructor"
            _LOGGER.debug(msg)
            brief= f"Instantiate {recipe.base} object using default constructor."
            ret_param = Param(
                "int32_t", "", "Object handle if successful and -1 for exception handling.")

        elif recipe.name == "del":
            # Default destructor
            recipe.what = "destructor"
            msg = f"   generating {func_name!r} -> default destructor"
            _LOGGER.debug(msg)
            args = [Param("int32_t", "handle", f"Handle to {recipe.base} object.")]
            brief= f"Delete {recipe.base} object."
            ret_param = Param(
                "int32_t", "", "Zero for success and -1 for exception handling.")

        else:
            msg = f"Unable to resolve recipe type for {recipe.name!r}"
            _LOGGER.critical(msg)
            exit(1)

        if recipe.brief:
            brief = recipe.brief
        uses = [self._doxygen_tags.cxx_member(uu) for uu in recipe.uses]
        return Func(ret_param.p_type, func_name, ArgList(args), brief, cxx_member,
                    ret_param.description, None, uses)

    def _handle_crosswalk(
            self, what: str, crosswalk: dict, derived: dict[str, str]) -> str:
        """Crosswalk for object handle."""
        cabinet = None
        classes = list(self._clib_bases) + list(derived.keys())
        for base in classes:
            ret_type = what.replace(f"<{base}>", "<T>")
            if ret_type in crosswalk:
                ret_type = crosswalk[ret_type]
                cabinet = base
                break
        if cabinet:
            # successful crosswalk with cabinet object
            return cabinet

        msg = f"Failed crosswalk for handle type {what!r} using {classes}."
        _LOGGER.critical(msg)
        sys.exit(1)

    def _ret_crosswalk(
            self, what: str, derived: list[str]) -> tuple[Param, list[Param]]:
        """Crosswalk for return type."""
        what = what.removeprefix("virtual ")
        ret_key = what.removeprefix("const ").removesuffix(" const").rstrip("&")
        if ret_key in self._config.ret_type_crosswalk:
            ret_type = self._config.ret_type_crosswalk[ret_key]
            if what.startswith("const "):
                ret_type = f"const {ret_type}"
            if ret_type == "void":
                returns = Param(
                    "int32_t", "", "Zero for success or -1 for exception handling.")
                return returns, []
            if ret_type.endswith("char*"):
                # string expressions require special handling
                returns = Param(
                    "int32_t", "", "Actual length of string including string-terminating "
                    "null byte, \\0, or -1 for exception handling.")
                ret_type = ret_type.removeprefix("const ")
                buffer = [
                    Param("int32_t", "bufLen", "Length of reserved array.", "in"),
                    Param(ret_type, "buf", "Returned string value.", "out")]
                return returns, buffer
            if ret_type.endswith("*"):
                # return type involves pointer to reserved buffer
                returns = Param(
                    "int32_t", "",
                    "Actual length of value array or -1 for exception handling.")
                ret_type = ret_type.removeprefix("const ")
                buffer = [
                    Param("int32_t", "bufLen", "Length of reserved array.", "in"),
                    Param(ret_type, "buf", "Returned array value.", "out")]
                return returns, buffer
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), []

            raise NotImplementedError(f"Crosswalk not implemented for {what!r}.")

        if "shared_ptr" in what:
            # check for crosswalk with object from includes
            handle = self._handle_crosswalk(
                what, self._config.ret_type_crosswalk, derived)
            returns = Param(
                "int32_t", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, []

        msg = f"Failed crosswalk for return type {what!r}."
        _LOGGER.critical(msg)
        sys.exit(1)


    def _prop_crosswalk(self, par_list: list[Param]) -> list[Param]:
        """Crosswalk for argument type."""
        if not par_list:
            return []
        params = []
        for par in par_list:
            what = par.p_type
            par_key = what.removeprefix("const ").removesuffix(" const").rstrip("&")
            if par_key in self._config.par_type_crosswalk:
                if "vector<" in par_key:
                    params.append(
                        Param("int32_t", f"{par.name}Len",
                              f"Length of vector reserved for {par.name}.", "in"))
                elif par_key.endswith("*"):
                    params.append(
                        Param("int32_t", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", "in"))
                ret_type = self._config.par_type_crosswalk[par_key]
                if what.startswith("const "):
                    ret_type = f"const {ret_type}"
                params.append(Param(ret_type, par.name, par.description, par.direction))
            elif "shared_ptr" in par_key:
                handle = self._handle_crosswalk(
                    par_key, self._config.par_type_crosswalk, {})
                par_key = par_key.replace(handle, "T")
                if "vector<" in par_key:
                    params.append(
                        Param("int32_t", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", "in"))
                    description = f"Memory holding {handle} objects. "
                    description += par.description
                    ret_type = self._config.par_type_crosswalk[par_key]
                    if what.startswith("const "):
                        ret_type = f"const {ret_type}"
                    params.append(Param(ret_type, par.name, description.strip()))
                else:
                    description = f"Integer handle to {handle} object. "
                    description += par.description
                    ret_type = self._config.par_type_crosswalk[par_key]
                    if what.startswith("const "):
                        ret_type = f"const {ret_type}"
                    params.append(
                        Param(ret_type, par.name, description.strip(), par.direction))
            else:
                msg = f"Failed crosswalk for argument type {par_key!r}."
                _LOGGER.critical(msg)
                sys.exit(1)
        return params
