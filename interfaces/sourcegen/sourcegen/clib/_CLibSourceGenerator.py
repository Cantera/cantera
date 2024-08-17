# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import logging
from typing import List


from ._Config import Config

from .._dataclasses import HeaderFile
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser


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
                print(f"{recipe.name}: {tag_info.signature()}")
