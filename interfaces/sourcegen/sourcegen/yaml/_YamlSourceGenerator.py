"""
Generator for YAML output.

Used for illustration purposes: the `CLibSourceGenerator` is used to preprocess YAML
header specifications, which yields the `HeaderFile` objects used by this source
generator.
"""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging

from jinja2 import Environment, BaseLoader

from .._dataclasses import HeaderFile, CFunc
from .._SourceGenerator import SourceGenerator


_LOGGER = logging.getLogger()

class YamlSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating YAML summaries of header specifications."""

    def __init__(self, out_dir: str, config: dict, templates: dict) -> None:
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = config
        self._templates = templates

    def _write_yaml(self, headers: HeaderFile) -> None:
        """Parse header file and generate YAML output."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)

        definition = loader.from_string(self._templates["yaml-definition"])
        declarations = []
        for c_func, recipe in zip(headers.funcs, headers.recipes):
            msg = f"   scaffolding {c_func.name!r} implementation"
            _LOGGER.debug(msg)
            implements = ""
            if isinstance(c_func.implements, CFunc):
                implements = c_func.implements.short_declaration()
            declarations.append(
                definition.render(c_func=c_func,
                                  returns=c_func.returns, implements=implements,
                                  relates=c_func.uses, what=recipe.what))

        filename = headers.output_name(suffix="3.yaml")
        template = loader.from_string(self._templates["yaml-file"])
        output = template.render(filename=filename.name, header_entries=declarations)

        out = Path(self._out_dir) / "yaml" / filename.name
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n")

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        """Generate output."""
        for headers in headers_files:
            msg = f"  parsing functions in {headers.path.name!r}"
            _LOGGER.info(msg)
            self._write_yaml(headers)
