"""Generator for Julia source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import sys
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from ..dataclasses import Func, HeaderFile
from ..generator import SourceGenerator

from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()
_LOADER = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the JuliaSourceGenerator class."""

    c_type_crosswalk: dict[str, str]  #: C type crosswalks


class JuliaSourceGenerator(SourceGenerator):
    """The SourceGenerator for scaffolding Julia files for the Julia interface"""

    def __init__(self, out_dir: str, config: dict, templates: dict) -> None:
        if not out_dir:
            _LOGGER.critical("Non-empty string identifying output path required.")
            sys.exit(1)
        self._out_dir = Path(out_dir) / "src" / "generated"

        # use the typed config
        self._config = Config(**config)
        self._templates = templates

    def _crosswalk(self, c_type: str) -> str:
        """Crosswalk of C/Julia types."""
        if c_type not in self._config.c_type_crosswalk:
            msg = (f"Unable to crosswalk C type {c_type!r}: add an entry to the "
                   "'c_type_crosswalk' section of sourcegen/julia/config.yaml.")
            _LOGGER.critical(msg)
            sys.exit(1)
        return self._config.c_type_crosswalk[c_type]

    def _scaffold_func(self, func: Func) -> str:
        """Render the `ccall` wrapper for a single CLib function."""
        arg_names = [par.name for par in func.arglist]
        arg_types = [self._crosswalk(par.p_type) for par in func.arglist]

        template = _LOADER.from_string(self._templates["julia-ccall-func"])
        return template.render(
            name=func.name,
            ret_type=self._crosswalk(func.ret_type),
            arg_names=", ".join(arg_names),
            arg_types="({},)".format(", ".join(arg_types)) if arg_types else "()",
            call_args="".join(f", {name}" for name in arg_names))

    def _write_file(self, file_name: str, template_name: str, **kwargs) -> None:
        _LOGGER.info(f"  writing {file_name!r}")
        t_file = Path(__file__).parent / template_name
        template = _LOADER.from_string(t_file.read_text(encoding="utf-8"))
        contents = template.render(file_name=file_name, **kwargs)

        self._out_dir.joinpath(file_name).write_text(contents, encoding="utf-8")

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        self._out_dir.mkdir(parents=True, exist_ok=True)
        # delete any existing files, so that bindings for a header that is no
        # longer generated cannot linger and shadow the current CLib
        for f in self._out_dir.iterdir():
            f.unlink()

        generated_files = []
        for header_file in headers_files:
            file_name = f"lib{header_file.path.stem}.jl"
            self._write_file(
                file_name, "template_bindings.jl.j2",
                header_file=f"cantera_clib/{header_file.path.stem}.h",
                functions=[self._scaffold_func(func) for func in header_file.funcs])
            generated_files.append(file_name)

        self._write_file("_manifest.jl", "template_manifest.jl.j2",
                         files=generated_files)
