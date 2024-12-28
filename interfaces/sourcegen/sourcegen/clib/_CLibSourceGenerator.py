"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging
# from typing import List, Tuple, Dict

# from jinja2 import Environment, BaseLoader

from ._Config import Config

from .._dataclasses import HeaderFile #, ArgList, Param, Recipe, Func
from .._SourceGenerator import SourceGenerator
# from .._TagFileParser import TagFileParser, TagDetails, tag_lookup


_logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    def __init__(self, out_dir: str, config: dict, templates: dict):
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config.from_parsed(**config)  # typed config
        self._templates = templates

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output"""
        raise NotImplementedError("Work in progress ...")
        # self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        # for header in headers_files:
        #     self._parse_header(header)
