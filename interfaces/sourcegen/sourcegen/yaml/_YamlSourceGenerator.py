"""
Generator for YAML output.

Used for debugging purposes only.
"""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging

from jinja2 import Environment, BaseLoader

from .._dataclasses import HeaderFile, CFunc
from .._SourceGenerator import SourceGenerator


_logger = logging.getLogger()

class YamlSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating YAML summaries of header specifications."""

    def __init__(self, out_dir: str, config: dict, templates: dict):
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = config
        self._templates = templates

    def _parse_header(self, header: HeaderFile):
        """Parse header file and generate output."""
        loader = Environment(loader=BaseLoader)

        filename = header.path.name
        template = loader.from_string(self._templates["yaml-preamble"])
        preamble = template.render(filename=filename)
        definition = loader.from_string(self._templates["yaml-definition"])
        declarations = []
        for c_func in header.funcs:

            implements = ""
            if isinstance(c_func.implements, CFunc):
                implements = c_func.implements.short_declaration()
            declarations.append(
                definition.render(c_func=c_func,
                                  returns=c_func.returns, implements=implements,
                                  relates=c_func.relates))

        guard = f"__{filename.upper().replace('.', '_')}__"
        template = loader.from_string(self._templates["yaml-file"])
        output = template.render(
            preamble=preamble, header_entries=declarations, guard=guard)

        if self._out_dir:
            out = Path(self._out_dir) / "yaml" / filename
            _logger.info(f"  writing {filename!r}")
            if not out.parent.exists():
                out.parent.mkdir(parents=True, exist_ok=True)
            with open(out, "wt", encoding="utf-8") as stream:
                stream.write(output)
                stream.write("\n")
        else:
            print(output)

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output."""
        for header in headers_files:
            _logger.info(f"  parsing functions in {header.path.name!r}:")
            self._parse_header(header)
