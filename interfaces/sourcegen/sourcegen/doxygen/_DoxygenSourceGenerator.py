# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from typing import List, Dict

from ._Config import Config

from .._dataclasses import HeaderFile, AnnotatedFunc
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser


class DoxygenSourceGenerator(SourceGenerator):
    """The SourceGenerator for referencing CLib functions to doxygen information."""

    def __init__(self, out_dir: str, config: dict):
        if out_dir:
            raise ValueError("Doxygen generator does not use output folder.")

        # use the typed config
        self._config = Config.from_parsed(config)
        self._doxygen = TagFileParser.from_parsed(config)

    def generate_source(self, headers_files: List[HeaderFile]):
        documented_funcs: Dict[str, AnnotatedFunc] = {}
        undocumented_funcs: List[str] = []

        for header_file in headers_files:
            doxy_funcs = list(map(self._doxygen.annotated_func, header_file.funcs))
            documented_funcs.update((f.name, f) for f in doxy_funcs if f is not None)
            undocumented_funcs.extend([_ for _ in doxy_funcs if _ is None])

        for func in documented_funcs.values():
            print(f"{func.name:<18}: {func.cxx_name:<20} "
                  f"{func.cxx_anchorfile} ({func.cxx_anchor})")
