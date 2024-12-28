"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging
# from typing import List, Tuple, Dict

from jinja2 import Environment, BaseLoader

from ._Config import Config

from .._dataclasses import HeaderFile #, ArgList, Param, Recipe, Func
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser #, TagDetails, tag_lookup


_logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    @staticmethod
    def _get_bases(headers_files: list[HeaderFile]) -> list[str]:
        bases = set()
        for headers in headers_files:
            for recipe in headers.recipes:
                bases |= set([recipe.base] + recipe.uses + recipe.parents)
        return list(bases)

    def _parse_header(self, header: HeaderFile):
        """Parse header file and generate output."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for recipe in header.recipes:
            _logger.info(f"    recipe for {recipe.name!r}")
        #     declaration, annotations, _ = self.build_declaration(recipe)
        #     declarations.append(
        #         template.render(declaration=declaration, annotations=annotations))

        filename = header.output_name(suffix=".h", auto="3")
        guard = f"__{filename.name.upper().replace('.', '_')}__"
        template = loader.from_string(self._templates["clib-header-file"])
        output = template.render(
            name=filename.stem, guard=guard, header_entries=declarations)

        if self._out_dir:
            out = Path(self._out_dir) / "include" / filename.name
            _logger.info(f"  writing {filename!r}")
            if not out.parent.exists():
                out.parent.mkdir(parents=True)
            with open(out, "wt", encoding="utf-8") as stream:
                stream.write(output)
                stream.write("\n")
        else:
            print(output)

    def __init__(self, out_dir: str, config: dict, templates: dict):
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config.from_parsed(**config)  # typed config
        self._templates = templates

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output."""
        self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        for header in headers_files:
            _logger.info(f"  parsing recipes in {header.path.name!r}:")
            self._parse_header(header)
