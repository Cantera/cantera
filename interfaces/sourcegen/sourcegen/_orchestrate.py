# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from distutils.command.config import LANG_EXT
import importlib
import inspect
from pathlib import Path
import ruamel.yaml

from ._dataclasses import HeaderFile
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
            config = ruamel.yaml.safe_load(config_file)

    ignore_files: list[str] = config.get("ignore_files", [])
    ignore_functions: dict[str, list[str]] = config.get("ignore_funcs", {})

    files = (HeaderFile.parse(f, ignore_functions.get(f.name, []))
        for f in _clib_path.glob("*.h") if f.name not in ignore_files)
    # removes instances where HeaderFile.parse() returned None
    files = list(filter(lambda p: p is not None and p.path != _clib_defs_path, files))

    # find and instantiate the language-specific SourceGenerator
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(Path(out_dir), config)

    scaffolder.generate_source(files)
