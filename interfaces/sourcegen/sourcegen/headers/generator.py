"""Resolve recipes based on specification file contents."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from .tagfiles import TagFileParser
from .._dataclasses import HeaderFile, Param, ArgList, CFunc, Recipe
from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the HeaderGenerator class"""

    ret_type_crosswalk: dict[str, str]  #: Return type cross-walks

    prop_type_crosswalk: dict[str, str]  #: Parameter type cross-walks


class HeaderGenerator:
    """Convert recipes to header file information."""

    _clib_bases: list[str] = None  #: list of bases provided via YAML configurations

    def __init__(self, config: dict, templates: dict, bases: list[str]) -> None:
        self._config = Config(**config)
        self._templates = templates
        self._doxygen_tags = None
        self._clib_bases = bases

    def resolve_tags(self, headers_files: list[HeaderFile]) -> None:
        """Resolve recipe information based on doxygen tags."""
        def get_bases() -> tuple[list[str], list[str]]:
            bases = set()
            classes = set()
            for headers in headers_files:
                bases.add(headers.base)
                for recipe in headers.recipes:
                    classes.update(recipe.bases)
            return sorted(bases), sorted(classes)

        self._clib_bases, classes = get_bases()
        self._doxygen_tags = TagFileParser(classes)

        for headers in headers_files:
            msg = f"  resolving recipes in {headers.path.name!r}:"
            _LOGGER.info(msg)
            c_funcs = []
            for recipe in headers.recipes:
                c_funcs.append(self.resolve_recipe(recipe))
            headers.funcs = c_funcs

    def resolve_recipe(self, recipe: Recipe) -> CFunc:
        """Build CLib header from recipe and doxygen annotations."""
        def merge_params(
                implements: str, cxx_member: CFunc | Param
            ) -> tuple[list[Param], CFunc]:
            """Create preliminary CLib argument list."""
            obj_handle = []
            if "::" in implements:
                # If class method, add handle as first parameter
                what = implements.split("::")[0]
                obj_handle.append(
                    Param("int", "handle", f"Handle to queried {what} object."))
            if isinstance(cxx_member, Param):
                if recipe.what.endswith("setter"):
                    return obj_handle + [cxx_member], cxx_member
                return obj_handle, cxx_member

            if "(" not in implements:
                return obj_handle + cxx_member.arglist.params, cxx_member

            # Signature may skip C++ default parameters
            args_short = CFunc.from_str(implements).arglist
            if len(args_short) < len(cxx_member.arglist):
                cxx_arglist = ArgList(cxx_member.arglist[:len(args_short)])
                cxx_member = CFunc(cxx_member.ret_type, cxx_member.name,
                                   cxx_arglist, cxx_member.brief, cxx_member.implements,
                                   cxx_member.returns, cxx_member.base, cxx_member.uses)

            return obj_handle + cxx_member.arglist.params, cxx_member

        func_name = f"{recipe.prefix}_{recipe.name}"
        reserved = ["cabinetSize", "parentHandle",
                    "getCanteraError", "setLogWriter", "setLogCallback",
                    "clearStorage", "resetStorage"]
        if recipe.name in reserved:
            recipe.what = "reserved"
            loader = Environment(loader=BaseLoader)
            msg = f"   generating {func_name!r} -> {recipe.what}"
            _LOGGER.debug(msg)
            header = loader.from_string(
                self._templates[f"clib-reserved-{recipe.name}-h"]
                ).render(base=recipe.base, prefix=recipe.prefix)
            return CFunc.from_str(header, brief=recipe.brief)

        # Ensure that all functions/methods referenced in recipe are detected correctly
        bases = recipe.bases
        if not recipe.implements:
            recipe.implements = self._doxygen_tags.detect(recipe.name, bases)
        elif recipe.base and "::" not in recipe.implements:
            parts = list(recipe.implements.partition("("))
            parts[0] = self._doxygen_tags.detect(parts[0], bases)
            recipe.implements = "".join(parts)
        recipe.uses = [self._doxygen_tags.detect(uu.split("(")[0], bases, False)
                       for uu in recipe.uses]

        cxx_member = None
        ret_param = Param("void")
        args = []
        brief = ""

        if recipe.implements:
            cxx_member = self._doxygen_tags.cxx_member(
                recipe.implements, recipe.what.endswith("setter"))

            if cxx_member.base in recipe.derived:
                # Use alternative prefix for class specialization
                recipe.prefix = recipe.derived.get(cxx_member.base, recipe.prefix)
                func_name = f"{recipe.prefix}_{recipe.name}"

            msg = f"   generating {func_name!r} -> {recipe.implements}"
            _LOGGER.debug(msg)

            if isinstance(cxx_member, CFunc):
                # Convert C++ return type to format suitable for crosswalk:
                # Incompatible return parameters are buffered and appended to back
                ret_param, buffer_params = self._ret_crosswalk(
                    cxx_member.ret_type, recipe.derived)
                par_list, cxx_member = merge_params(recipe.implements, cxx_member)
                prop_params = self._prop_crosswalk(par_list)
                args = prop_params + buffer_params
                brief = cxx_member.brief
            elif recipe.what == "variable-setter":
                ret_param = Param("int")
                par_list, cxx_member = merge_params(recipe.implements, cxx_member)
                args = self._prop_crosswalk(par_list)
                brief = cxx_member.description
            else:
                # Variable getter
                prop_params, cxx_member = merge_params(recipe.implements, cxx_member)
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
            elif any(recipe.implements.startswith(base) for base in recipe.bases):
                recipe.what = "method"
            else:
                msg = f"Unable to auto-detect function type for recipe {recipe.name!r}."
                _LOGGER.critical(msg)
                exit(1)

        elif recipe.what == "noop":
            # No operation
            msg = f"   generating {func_name!r} -> no-operation"
            _LOGGER.debug(msg)
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            brief = "No operation."
            ret_param = Param("int", "", "Always zero.")

        elif recipe.name == "new":
            # Default constructor
            recipe.what = "constructor"
            msg = f"   generating {func_name!r} -> default constructor"
            _LOGGER.debug(msg)
            brief= f"Instantiate {recipe.base} object using default constructor."
            ret_param = Param(
                "int", "", "Object handle if successful and -1 for exception handling.")

        elif recipe.name == "del":
            # Default destructor
            recipe.what = "destructor"
            msg = f"   generating {func_name!r} -> default destructor"
            _LOGGER.debug(msg)
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            brief= f"Delete {recipe.base} object."
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")

        elif recipe.code:
            # Custom code
            return CFunc.from_str(recipe.code, brief=recipe.brief)

        else:
            msg = f"Unable to resolve recipe type for {recipe.name!r}"
            _LOGGER.critical(msg)
            exit(1)

        if recipe.brief:
            brief = recipe.brief
        uses = [self._doxygen_tags.cxx_member(uu) for uu in recipe.uses]
        return CFunc(ret_param.p_type, func_name, ArgList(args), brief, cxx_member,
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
        what = what.replace("virtual ", "")
        if what in self._config.ret_type_crosswalk:
            ret_type = self._config.ret_type_crosswalk[what]
            if ret_type == "void":
                returns = Param(
                    "int", "", "Zero for success or -1 for exception handling.")
                return returns, []
            if ret_type == "char*":
                # string expressions require special handling
                returns = Param(
                    "int", "", "Actual length of string including string-terminating "
                    "null byte, \\0, or -1 for exception handling.")
                buffer = [
                    Param("int", "bufLen", "Length of reserved array.", "in"),
                    Param(ret_type, "buf", "Returned string value.", "out")]
                return returns, buffer
            if ret_type.endswith("*"):
                # return type involves pointer to reserved buffer
                returns = Param(
                    "int", "",
                    "Actual length of value array or -1 for exception handling.")
                buffer = [
                    Param("int", "bufLen", "Length of reserved array.", "in"),
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
                "int", "",
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
            if what in self._config.prop_type_crosswalk:
                if "vector<" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of vector reserved for {par.name}.", "in"))
                elif what.endswith("* const") or what.endswith("double*"):
                    direction = "in" if what.startswith("const") else "out"
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", direction))
                ret_type = self._config.prop_type_crosswalk[what]
                params.append(Param(ret_type, par.name, par.description, par.direction))
            elif "shared_ptr" in what:
                handle = self._handle_crosswalk(
                    what, self._config.prop_type_crosswalk, {})
                if "vector<" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", "in"))
                    description = f"Memory holding {handle} objects. "
                    description += par.description
                    params.append(Param("const int*", par.name, description.strip()))
                else:
                    description = f"Integer handle to {handle} object. "
                    description += par.description
                    params.append(
                        Param("int", par.name, description.strip(), par.direction))
            else:
                msg = f"Failed crosswalk for argument type {what!r}."
                _LOGGER.critical(msg)
                sys.exit(1)
        return params
