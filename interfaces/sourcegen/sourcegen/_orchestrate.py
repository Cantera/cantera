# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib
from pathlib import Path
import re
import ruamel.yaml

from ._types import Func


def _call_generate_source(generator, incl_file, ignore):
    print('  ' + incl_file.name)

    ct = incl_file.read_text()

    matches = re.finditer(r'CANTERA_CAPI.*?;', ct, re.DOTALL)
    c_functions = [re.sub(r'\s+', ' ', m.group()) for m in matches]

    if not c_functions:
        return

    parsed = map(Func.parse, c_functions)

    if ignore:
        print(f'    ignoring ' + str(ignore))

    parsed = [f for f in parsed if f.name not in ignore]

    generator.generate_source(incl_file, parsed)


def generate_source(lang: str, out_dir: str):
    print('Generating source files...')

    my_path = Path(__file__).parent
    clib_path = my_path.joinpath('../../../include/cantera/clib').resolve()
    config_path = my_path.joinpath(lang, 'config.yaml').resolve()

    with config_path.open() as config_file:
        config = ruamel.yaml.safe_load(config_file)

    module = importlib.import_module(__package__  + '.' + lang)
    generator = module.SourceGenerator(Path(out_dir), config)

    ignore = config['ignore']

    for incl_file in clib_path.glob('*.h'):
        if incl_file.name not in ignore or ignore[incl_file.name]:
            _call_generate_source(generator, incl_file, ignore.get(incl_file.name, []))

    generator.finalize()
