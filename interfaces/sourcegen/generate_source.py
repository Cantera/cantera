import importlib
import os
import re
import sys

import ruamel.yaml

my_path = sys.path[0] + '/'


def split_param(param_string):
    parts = param_string.strip().rsplit(' ', 1)
    if len(parts) == 2:
        return tuple(parts)


def parse(c_function):
    lparen = c_function.index('(')
    rparen = c_function.index(')')
    front = c_function[0:lparen].split()
    
    params = (split_param(p) for p in c_function[lparen+1:rparen].split(','))
    params = [s for s in params if s]
    
    ret_type = front[-2]
    name = front[-1]
    return ret_type, name, params


def call_generate_source(incl_file, ignore):
    print('  ' + incl_file.path)

    with open(incl_file.path, 'r') as f:
        ct = f.read()

    matches = re.finditer(r'CANTERA_CAPI.*?;', ct, re.DOTALL)
    c_functions = [re.sub(r'\s+', ' ', m.group()) for m in matches]

    if not c_functions:
        return

    parsed = (parse(f) for f in c_functions)

    if ignore:
        print(f'    ignoring ' + str(ignore))

    parsed = [(ret_type, name, params) for ret_type, name, params in parsed if name not in ignore]

    generator.generate_source(incl_file, parsed)


lang, out_dir = sys.argv[1:]

print('Gnenerating source files...')

with open(my_path + lang + '/config.yaml', 'r') as config_file:
   config = ruamel.yaml.safe_load(config_file)

generator = importlib.import_module(lang).SourceGenerator(out_dir, config)

ignore = config['ignore']

for incl_file in os.scandir(my_path + '../../include/cantera/clib'):
    if incl_file.name not in ignore or ignore[incl_file.name]:
        call_generate_source(incl_file, ignore.get(incl_file.name, []))

generator.finalize()