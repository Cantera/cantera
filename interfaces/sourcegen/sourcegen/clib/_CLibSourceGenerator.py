"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import logging
from typing import List, Tuple, Dict


from ._Config import Config

from .._dataclasses import HeaderFile, ArgList, Param, Recipe, Func
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser, TagDetails, tag_lookup


logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    @staticmethod
    def _get_bases(headers_files: List[HeaderFile]) -> List[str]:
        bases = set()
        for headers in headers_files:
            for recipe in headers.recipes:
                bases |= set([recipe.base] + recipe.uses)
        return list(bases)

    @staticmethod
    def _build_annotation(details: TagDetails, ret: Param,
                          params: List[Param], relates="") -> str:
        """Build annotation block."""
        msg = ["/**", details.briefdescription, ""]

        def param(item: Param):
            ret = "@param"
            if item.direction:
                ret += f"[{item.direction}]"
            ret += f" {item.name}"
            if item.description:
                return f"{ret:<20} {item.description}"
            return f"{ret:<20} Undocumented."

        msg += [param(_) for _ in params]
        if ret.description:
            msg += [f"{'@returns':<20} {ret.description}"]
        arglist = ArgList.from_xml(details.arglist)
        msg += ["", f"@implements {details.qualified_name}{arglist.short_str()}"]
        if relates:
            msg += [f"@relates {relates}"]
        return "\n * ".join(msg).strip() + "\n */"

    def _handle_crosswalk(self, what: str, crosswalk: Dict) -> str:
        """Crosswalk for object handle."""
        cabinet = None
        for base in self._config.cabinets:
            ret_type = what.replace(f"<{base}>", "<T>")
            if ret_type in crosswalk:
                ret_type = crosswalk[ret_type]
                cabinet = base
                break
        if cabinet:
            # successful crosswalk with cabinet object
            return cabinet

        logging.critical("Failed crosswalk for handle type '%s'.", what)
        sys.exit(1)

    def _ret_crosswalk(self, what: str) -> Tuple[Param, List[Param], List[str]]:
        """Crosswalk for return type."""
        what = what.replace("virtual ", "")
        if what in self._config.ret_type_crosswalk:
            ret_type = self._config.ret_type_crosswalk[what]
            if ret_type == "char*":
                # string expressions require special handling
                returns = Param(
                    "int", "",
                    "Actual length of string or -1 for exception handling.")
                buffer = [
                    Param("int", "lenBuf", "Length of reserved array.", "in"),
                    Param(ret_type, "charBuf", "Returned string value.", "out")]
                return returns, buffer, []
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), [], []

            # all other types require reserved buffers
            what = "vector" if what.startswith("vector") else "array"
            returns = Param(
                "int", "", f"Actual length of {what} or -1 for exception handling.")
            buffer = [
                Param("int", "lenBuf", "Length of reserved array.", "in"),
                Param(ret_type, "valueBuf", f"Returned {what} value.", "out")]
            return returns, buffer, []

        if what.startswith("shared_ptr"):
            # check for crosswalk with object from cabinets
            handle = self._handle_crosswalk(what, self._config.ret_type_crosswalk)
            returns = Param(
                "int", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, [], [handle]

        logging.critical("Failed crosswalk for return type '%s'.", what)
        sys.exit(1)

    def _prop_crosswalk(self, par_list: List[Param]) -> Tuple[List[Param], List[str]]:
        """Crosswalk for argument type."""
        if not par_list:
            return [], []
        params = []
        cabinets = set()
        for par in par_list:
            what = par.p_type
            if what in self._config.prop_type_crosswalk:
                if "vector" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of vector reserved for {par.name}.", "in"))
                ret_type = self._config.prop_type_crosswalk[what]
                params.append(Param(ret_type, par.name, par.description, par.direction))
            elif "shared_ptr" in what:
                handle = self._handle_crosswalk(what, self._config.prop_type_crosswalk)
                cabinets |= {handle}
                params.append(Param("int", par.name, par.description, par.direction))
            else:
                logging.critical("Failed crosswalk for argument type '%s'.", what)
                sys.exit(1)
        return params, list[cabinets]

    def build_declaration(self, recipe: Recipe) -> Tuple[str, str, List[str]]:
        """Build strings containing declaration and annotation."""
        def merge_params(implements, details: TagDetails) -> List[Param]:
            # If class method, add handle as first parameter
            args_merged = []
            if "::" in implements:
                what = implements.split("::")[0]
                args_merged.append(
                    Param("int", "handle", f"Handle to {what} object."))
            # Merge parameters from signature, doxygen info and doxygen details
            args_used = ArgList.from_xml(details.arglist).params  # from doxygen
            if "(" in implements:
                args_short = Func.from_str(implements).params  # from recipe
                args_used = args_used[:len(args_short)]
            args_annotated = details.parameterlist  # from documentation
            for arg in args_used:
                for desc in args_annotated:
                    if arg.name == desc.name:
                        args_merged.append(
                            Param(arg.p_type, arg.name,
                                  desc.description, desc.direction, arg.default))
                        break
                else:
                    args_merged.append(arg)
            return args_merged

        if recipe.implements:
            tag_info = self._doxygen_tags.tag_info(recipe.implements)
            details = tag_lookup(tag_info)

            # convert XML return type to format suitable for crosswalk
            ret_type = Param.from_xml(details.type).p_type
            ret_param, buffer_params, cabinets = self._ret_crosswalk(ret_type)
            par_list = merge_params(recipe.implements, details)
            prop_params, prop_cabinets = self._prop_crosswalk(par_list)
            cabinets += prop_cabinets
            all_params = prop_params + buffer_params
            c_func = Func("", ret_param.p_type, recipe.name, all_params, "")
            declaration = f"{c_func.declaration()};"
            annotations = self._build_annotation(
                details, ret_param, all_params, recipe.relates)

        elif recipe.what == "destructor":
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            details = TagDetails(
                "", "", "", "", "", "", f"Delete {recipe.base} handle.", "", args)
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")
            annotations = f"//! {details.briefdescription}"
            buffer_params = []
            cabinets = [recipe.base]
            declaration = f"{ret_param.p_type} {recipe.name}(...);"

        else:
            logger.critical("Unable to build declaration for '%s' with type '%s'.",
                            recipe.name, recipe.what)
            sys.exit(1)


        return declaration, annotations, cabinets

    def _parse_header(self, header: HeaderFile):
        for recipe in header.recipes:
            declaration, annotations, _ = self.build_declaration(recipe)
            print(f"{annotations}\n{declaration}\n")

    def __init__(self, out_dir: str, config: dict):
        self._out_dir = out_dir or None
        self._config = Config.from_parsed(**config)  # typed config

    def generate_source(self, headers_files: List[HeaderFile]):
        """Generate output"""
        self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        for header in headers_files:
            self._parse_header(header)
