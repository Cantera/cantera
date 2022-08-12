# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
import inspect
from pathlib import Path
import ruamel.yaml

from ._dataclasses import HeaderFile
from ._SourceGenerator import SourceGenerator


def generate_source(lang: str, out_dir: str):
    print("Generating source files...")

    my_path = Path(__file__).parent
    clib_path = my_path.joinpath("../../../include/cantera/clib").resolve()
    config_path = my_path.joinpath(lang, "config.yaml").resolve()

    with config_path.open() as config_file:
        config = ruamel.yaml.safe_load(config_file)

    # find and instantiate the language-specific SourceGenerator
    module = importlib.import_module(__package__ + "." + lang)
    _, scaffolder_type = inspect.getmembers(module,
        lambda m: inspect.isclass(m) and issubclass(m, SourceGenerator))[0]
    scaffolder: SourceGenerator = scaffolder_type(Path(out_dir), config)

    ignore_files: list[str] = config.get("ignore_files", [])
    ignore_functions: dict[str, list[str]] = config.get("ignore_functions", {})

    files = (HeaderFile.parse(f, ignore_functions.get(f.name, []))
        for f in clib_path.glob("*.h") if f.name not in ignore_files)
    # removes instances where HeaderFile.parse() returned None
    files = list(filter(None, files))

    scaffolder.generate_source(files)
