import importlib
import os
import os.path
import re

import ruamel.yaml

from ._types import *


def _call_generate_source(generator, incl_file, ignore):
    print('  ' + incl_file.path)

    with open(incl_file.path, 'r') as f:
        ct = f.read()

    matches = re.finditer(r'CANTERA_CAPI.*?;', ct, re.DOTALL)
    c_functions = [re.sub(r'\s+', ' ', m.group()) for m in matches]

    if not c_functions:
        return

    parsed = (Func.parse(f) for f in c_functions)

    if ignore:
        print(f'    ignoring ' + str(ignore))

    parsed = [f for f in parsed if f.name not in ignore]

    generator.generate_source(incl_file, parsed)


def generate_source(lang: str, out_dir: str):
    print('Gnenerating source files...')

    my_path = os.path.dirname(__file__) + '/'

    with open(my_path + lang + '/config.yaml', 'r') as config_file:
        config = ruamel.yaml.safe_load(config_file)

    generator = importlib.import_module(__package__  + '.' + lang).SourceGenerator(out_dir, config)

    ignore = config['ignore']

    for incl_file in os.scandir(os.path.normpath(my_path + '../../../include/cantera/clib')):
        if incl_file.name not in ignore or ignore[incl_file.name]:
            _call_generate_source(generator, incl_file, ignore.get(incl_file.name, []))

    generator.finalize()
