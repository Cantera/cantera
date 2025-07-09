# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
import logging
import sys
from typing import List, Dict

from ._HeaderFileParser import HeaderFileParser
from ._SourceGenerator import SourceGenerator
from ._helpers import read_config


_logger = logging.getLogger()

class CustomFormatter(logging.Formatter):
    """Minimalistic logging output"""

    def format(self, record):
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        return formatter.format(record)


def generate_source(lang: str, out_dir: str=""):
    """Main entry point of sourcegen."""
    loghandler = logging.StreamHandler(sys.stdout)
    loghandler.setFormatter(CustomFormatter())
    _logger.handlers.clear()
    _logger.addHandler(loghandler)
    _logger.setLevel(logging.DEBUG)
    _logger.info(f"Generating {lang!r} source files...")

    module = importlib.import_module(__package__ + "." + lang)
    root = Path(module.__file__).parent
    config = read_config(root / "config.yaml")
    templates = read_config(root / "templates.yaml")
    ignore_files: List[str] = config.pop("ignore_files", [])
    ignore_funcs: Dict[str, List[str]] = config.pop("ignore_funcs", {})

    files = HeaderFileParser.from_headers(ignore_files, ignore_funcs)

    # find and instantiate the language-specific SourceGenerator
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(out_dir, config, templates)

    scaffolder.generate_source(files)
    _logger.info("Done.")
