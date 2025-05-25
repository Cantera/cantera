"""Source generator for creating Cantera interface code."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
import logging
import sys
import argparse

from .headers import HeaderFileParser
from .generator import SourceGenerator
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

    if not lang:
        _LOGGER.critical("Aborting: sourcegen requires API language specification. "
                         "Use option --help for more information.")
        exit(1)

    if not out_dir:
        _LOGGER.critical("Aborting: sourcegen requires OUTPUT folder information. "
                         "Use option --help for more information.")
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


def create_argparser():
    """Sourcegen argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Source generator for creating Cantera interface code."),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="show additional logging output")
    parser.add_argument(
        "--api", choices=["clib", "csharp", "yaml"],
        help="language of generated Cantera API code")
    parser.add_argument(
        "--output", default="",
        help="specifies the OUTPUT folder name")

    return parser


def main():
    """Sourcegen CLI."""
    parser = create_argparser()
    args = parser.parse_args()
    lang = args.api
    output = args.output
    verbose = args.verbose
    generate_source(lang, output, verbose=verbose)
