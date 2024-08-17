"""Generator for annotated YAML summaries of CLib header files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging
from typing import List, Dict, Union
try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml

from ._Config import Config
from ._dataclasses import AnnotatedFunc

from .._dataclasses import HeaderFile, Func
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser
from .._HeaderFileParser import doxygen_func


logger = logging.getLogger()

class YamlSourceGenerator(SourceGenerator):
    """The SourceGenerator for referencing CLib functions to doxygen information."""

    def annotated_func(self, parsed: Func) -> Union[AnnotatedFunc, None]:
        """Match function with doxygen tag information."""
        comments = parsed.annotations

        implements = doxygen_func("@implements", comments)
        if not implements:
            return None

        tag_info = self._doxygen_tags.tag_info(implements)
        if not tag_info:
            logging.error("Unable to retrieve tag info for '%s'", implements)
            return None

        return AnnotatedFunc(*parsed,
                             implements,
                             doxygen_func("@relates", comments),
                             *tag_info)

    def __init__(self, out_dir: str, config: dict):
        self._out_dir = out_dir or None

        # use the typed config
        self._config = Config.from_parsed(**config)
        self._doxygen_tags = TagFileParser(self._config.bases)

    def generate_source(self, headers_files: List[HeaderFile]):
        """Generate output."""
        annotated_map: Dict[str, List[Dict[str, str]]] = {}

        for header_file in headers_files:
            name = header_file.path.name
            annotated_map[name] = []
            for func in list(map(self.annotated_func, header_file.funcs)):
                if func is not None:
                    annotated_map[name].append(func)

        emitter = yaml.YAML()
        emitter.width = 80
        emitter.register_class(AnnotatedFunc)
        if self._out_dir:
            out = Path(self._out_dir) / "interop.yaml"
            logger.info("  writing %s", out.name)
            with open(out, "wt", encoding="utf-8") as stream:
                emitter.dump(annotated_map, stream)
        else:
            emitter.dump(annotated_map, sys.stdout)
