"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import logging
from typing import List, Tuple, Dict


from ._Config import Config

from .._dataclasses import HeaderFile, ArgList, Param, Recipe
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser, TagDetails, tag_lookup


logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    # _doxygen_tags:

    @staticmethod
    def _get_bases(headers_files: List[HeaderFile]) -> List[str]:
        bases = set()
        for headers in headers_files:
            for recipe in headers.recipes:
                bases |= set([recipe.base] + recipe.uses)
        return list(bases)

    @staticmethod
    def _build_annotation(details: TagDetails, relates="") -> str:
        """Build annotation block."""
        msg = ["", details.briefdescription]

        def param(item: Param):
            ret = "@param"
            if item.direction:
                ret += f"[{item.direction}]"
            ret += f" {item.name}"
            return f"{ret:<20} {item.description}"

        msg += [param(_) for _ in details.parameterlist]
        arglist = ArgList.from_xml(details.arglist)
        msg += [f"@implements {details.qualified_name}{arglist.short_str()}"]
        if relates:
            msg += [f"@relates {relates}"]
        return "\n//! ".join(msg).strip()

    def _handle_crosswalk(self, what) -> str:
        """Crosswalk for object handle."""
        cabinet = None
        for base in self._config.cabinets:
            ret_type = what.replace(f"<{base}>", "<T>")
            if ret_type in self._config.ret_type_crosswalk:
                ret_type = self._config.ret_type_crosswalk[ret_type]
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
                    "int", "", "", "",
                    "Actual length of string or -1 for exception handling.")
                buffer = [
                    Param("int", "lenBuf", "", "in", "Length of reserved array."),
                    Param(ret_type, "charBuf", "", "out", "Returned string value.")]
                return returns, buffer, []
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), [], []

            # all other types require reserved buffers
            what = "vector" if what.startswith("vector") else "array"
            returns = Param(
                "int", "", "", "",
                f"Actual length of {what} or -1 for exception handling.")
            buffer = [
                Param("int", "lenBuf", "", "in", "Length of reserved array."),
                Param(ret_type, "valueBuf", "", "out", f"Returned {what} value.")]
            return returns, buffer, []

        if what.startswith("shared_ptr"):
            # check for crosswalk with object from cabinets
            handle = self._handle_crosswalk(what)
            returns = Param(
                "int", "", "", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, [], [handle]

        logging.critical("Failed crosswalk for return type '%s'.", what)
        sys.exit(1)

    def build_declaration(self, recipe: Recipe) -> Tuple[str, str, List[str]]:
        """Build strings containing declaration and annotation."""
        if recipe.implements:
            tag_info = self._doxygen_tags.tag_info(recipe.implements)
            details = tag_lookup(tag_info)

            # convert XML return type to format suitable for crosswalk
            ret_type = Param.from_xml(details.type).p_type
            ret_param, buffer_params, cabinets = self._ret_crosswalk(ret_type)
            annotations = self._build_annotation(details, recipe.relates)

        elif recipe.what == "destructor":
            details = TagDetails(
                "", "", "", "", "", "", f"Delete {recipe.base} handle.", "",
                [Param("int", "handle", "", "", f"Handle to {recipe.base} object.")])
            ret_param = Param(
                "int", "", "", "",
                "Zero for success and -1 for exception handling.")
            annotations = f"//! {details.briefdescription}"
            buffer_params = []
            cabinets = [recipe.base]

        else:
            logger.critical("Unable to build declaration for '%s'.", recipe.name)
            sys.exit(1)

        declaration = f"{ret_param.p_type} {recipe.name}(...);"

        return declaration, annotations, cabinets

    def _parse_header(self, header: HeaderFile):
        for recipe in header.recipes:
            declaration, annotations, _ = self.build_declaration(recipe)
            print(f"\n{annotations}\n{declaration}")

    def __init__(self, out_dir: str, config: dict):
        self._out_dir = out_dir or None
        self._config = Config.from_parsed(**config)  # typed config

    def generate_source(self, headers_files: List[HeaderFile]):
        """Generate output"""
        self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        for header in headers_files:
            self._parse_header(header)
