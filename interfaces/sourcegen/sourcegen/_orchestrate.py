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


logger = logging.getLogger()

class CustomFormatter(logging.Formatter):
    """Minimalistic logging output"""

    def format(self, record):
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        return formatter.format(record)


def generate_source(lang: str, out_dir: str=""):
    """Main entry point of sourcegen."""
    loghandler = logging.StreamHandler(sys.stdout)
    loghandler.setFormatter(CustomFormatter())
    logger.handlers.clear()
    logger.addHandler(loghandler)
    logger.setLevel(logging.DEBUG)
    logger.info(f"Generating {lang!r} source files...")

    module = importlib.import_module(__package__ + "." + lang)
    root = Path(module.__file__).parent
    config = read_config(root.joinpath("config.yaml"))
    templates = read_config(root.joinpath("templates.yaml"))
    ignore_files: List[str] = config.pop("ignore_files", [])
    ignore_funcs: Dict[str, List[str]] = config.pop("ignore_funcs", {})

    if lang == 'clib':
        files = HeaderFileParser.parse_yaml(ignore_files, ignore_funcs)
    else:
        files = HeaderFileParser.parse_h(ignore_files, ignore_funcs)

    # find and instantiate the language-specific SourceGenerator
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(out_dir, config, templates)

    scaffolder.generate_source(files)
    logger.info("Done.")
