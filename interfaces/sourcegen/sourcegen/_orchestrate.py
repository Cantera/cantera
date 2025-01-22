# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
import logging
import sys

from ._HeaderFileParser import HeaderFileParser
from ._SourceGenerator import SourceGenerator
from .clib import CLibSourceGenerator
from ._helpers import read_config


_LOGGER = logging.getLogger()

class CustomFormatter(logging.Formatter):
    """Minimalistic logging output"""

    def format(self, record):
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        return formatter.format(record)


def generate_source(lang: str, out_dir: str, verbose: bool = False) -> None:
    """Main entry point of sourcegen."""
    loghandler = logging.StreamHandler(sys.stdout)
    loghandler.setFormatter(CustomFormatter())
    _LOGGER.handlers.clear()
    _LOGGER.addHandler(loghandler)
    _LOGGER.setLevel(logging.DEBUG if verbose else logging.INFO)

    if not out_dir:
        _LOGGER.critical("Aborting: sourcegen requires output folder information.")
        exit(1)

    module = importlib.import_module(__package__ + "." + lang)
    root = Path(module.__file__).parent
    config = read_config(root / "config.yaml")
    templates = read_config(root / "templates.yaml")
    ignore_files: list[str] = config.pop("ignore_files", [])
    ignore_funcs: dict[str, list[str]] = config.pop("ignore_funcs", {})

    msg = f"Starting sourcegen for {lang!r} API"
    _LOGGER.info(msg)

    if lang == "clib":
        # prepare for generation of CLib headers in main processing step
        files = HeaderFileParser.headers_from_yaml(ignore_files, ignore_funcs)
    elif lang == "csharp":
        # csharp parses existing (traditional) CLib header files
        files = HeaderFileParser.headers_from_h(ignore_files, ignore_funcs)
    else:
        # generate CLib headers from YAML specifications as a preprocessing step
        files = HeaderFileParser.headers_from_yaml(ignore_files, ignore_funcs)
        clib_root = Path(__file__).parent / "clib"
        clib_config = read_config(clib_root / "config.yaml")
        clib_templates = read_config(clib_root / "templates.yaml")
        for key in ["ignore_files", "ignore_funcs"]:
            clib_config.pop(key)
        clib_scaffolder = CLibSourceGenerator(None, clib_config, clib_templates)
        clib_scaffolder.resolve_tags(files)

    # find and instantiate the language-specific SourceGenerator
    msg = f"Generating {lang!r} source files..."
    _LOGGER.info(msg)
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(out_dir, config, templates)

    scaffolder.generate_source(files)
    _LOGGER.info("Done.")
