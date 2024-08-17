"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import logging
from typing import List


from ._Config import Config

from .._dataclasses import HeaderFile, ArgList, Param
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
                bases |= set(recipe.bases)
        return list(bases)

    @staticmethod
    def _make_annotation(details: TagDetails, relates="") -> str:
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

    def __init__(self, out_dir: str, config: dict):
        self._out_dir = out_dir or None
        self._config = Config.from_parsed(**config)  # typed config

    def generate_source(self, headers_files: List[HeaderFile]):
        """Generate output"""
        self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        for headers in headers_files:
            for recipe in headers.recipes:
                if not recipe.implements:
                    print(f"{recipe.name}: N/A")
                    continue
                tag_info = self._doxygen_tags.tag_info(recipe.implements)
                tag_details = tag_lookup(tag_info)
                print(f"{recipe.name}:\n"
                      f"{self._make_annotation(tag_details, recipe.relates)}\n")
