# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
import logging
import sys

from ._HeaderFileParser import HeaderFileParser
from ._SourceGenerator import SourceGenerator
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

    # generate CLib headers from YAML specifications
    files = HeaderFileParser.headers_from_yaml(ignore_files, ignore_funcs)

    # find and instantiate the language-specific SourceGenerator
    msg = f"Generating {lang!r} source files..."
    _LOGGER.info(msg)
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(out_dir, config, templates)

    scaffolder.generate_source(files)
    _LOGGER.info("Done.")
