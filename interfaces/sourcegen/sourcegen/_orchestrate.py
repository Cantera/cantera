# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
from typing import List, Dict
import ruamel.yaml

from ._HeaderFileParser import HeaderFileParser
from ._SourceGenerator import SourceGenerator


_clib_path = Path(__file__).parent.joinpath("../../../include/cantera/clib").resolve()
_clib_defs_path = _clib_path.joinpath("clib_defs.h")

def generate_source(lang: str, out_dir: str):
    print("Generating source files...")

    module = importlib.import_module(__package__ + "." + lang)
    config_path = Path(module.__file__).parent.joinpath("config.yaml")

    config = {}
    if config_path.exists():
        with config_path.open() as config_file:
            reader = ruamel.yaml.YAML(typ="safe")
            config = reader.load(config_file)

    ignore_files: List[str] = config.get("ignore_files", [])
    ignore_funcs: Dict[str, List[str]] = config.get("ignore_funcs", {})

    files = (HeaderFileParser(f, ignore_funcs.get(f.name, [])).parse()
        for f in _clib_path.glob("*.h")
        if f != _clib_defs_path and f.name not in ignore_files)
    # removes instances where HeaderFile.parse() returned None
    files = list(filter(None, files))

    # find and instantiate the language-specific SourceGenerator
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(Path(out_dir), config)

    scaffolder.generate_source(files)
