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

from ..dataclasses import HeaderFile, Func, Param
from ..generator import SourceGenerator


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
            wraps = ""
            if isinstance(c_func.wraps, Func):
                wraps = c_func.wraps.short_declaration()
            elif isinstance(c_func.wraps, Param):
                wraps = c_func.wraps.long_str()
            uses = []
            for fcn in c_func.uses:
                uses.append(fcn.short_declaration())
            declarations.append(
                definition.render(c_func=c_func,
                                  returns=c_func.returns, wraps=wraps,
                                  uses=uses, what=recipe.what))

        filename = headers.output_name(suffix=".yaml")
        t_file = Path(__file__).parent / "template_output.yaml.j2"
        template = loader.from_string(t_file.read_text(encoding="utf-8"))
        output = template.render(filename=filename.name, header_entries=declarations)

        out = Path(self._out_dir) / filename.name
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n", encoding="utf-8")

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        """Generate output."""
        for headers in headers_files:
            msg = f"  parsing functions in {headers.path.name!r}"
            _LOGGER.info(msg)
            self._write_yaml(headers)
