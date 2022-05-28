from operator import indexOf
import os
import re
import sys

import ruamel.yaml

print('Gnenerating interop types...')

prolog = '    [DllImport(LibFile)]\n    public static extern'

type_map = {
    'const char*': 'string',
    'size_t': 'nuint',
    'char*': 'byte*'
}


def join_params(params):
    return ', '.join((param_name + ' ' + param_type for (param_name, param_type) in params))


def get_function_text(function):
    (ref_type, name, params, _) = function
    is_unsafe = any((param_type.endswith('*') for (param_type, _) in params))
    if is_unsafe:
        return f'{prolog} unsafe {ref_type} {name}({join_params(params)});'
    else:
        return f'{prolog} {ref_type} {name}({join_params(params)});'


def get_base_handle_text(handle):
    (name, del_clazz) = handle

    handle = f'''class {del_clazz} : CanteraHandle
{{
    protected override bool ReleaseHandle() =>
        Convert.ToBoolean(LibCantera.{name}(Value));
}}'''

    return handle


def get_derived_handle_text(derived):
    (derived, base) = derived

    derived = f'''class {derived} : {base} {{ }}'''

    return derived


def split_param(param_string):
    parts = param_string.strip().rsplit(' ', 1)
    if len(parts) == 2:
        return tuple(parts)


def parse(c_function):
    lparen = c_function.index('(')
    rparen = c_function.index(')')
    front = (c_function[0:lparen]).split()
    
    params = (split_param(p) for p in c_function[lparen+1:rparen].split(','))
    params = [s for s in params if s]
    
    ret_type = front[-2]
    name = front[-1]
    return (ret_type, name, params)


def convert(parsed):
    (ret_type, name, params) = parsed
    (clazz, method) = tuple(name.split('_', 1))
    
    del_clazz = None

    if clazz != 'ct':
        handle_clazz = clazz.capitalize() + 'Handle'
        
        # It’s not a “global” function, therefore:
        #   * It wraps a constructor and returns a handle, or
        #   * It wraps an instance method and takes the handle as the first parameter.
        if method.startswith('del'):
            del_clazz = handle_clazz
        elif method.startswith('new'):
            ret_type = handle_clazz
        else:
            (_, param_name) = params[0]
            params[0] = (handle_clazz, param_name)
            
    for (c_type, cs_type) in type_map.items():
        if (ret_type == c_type):
            ret_type = cs_type
            break
            
    for i in range(0, len(params)):
        (param_type, param_name) = params[i]
        
        for (c_type, cs_type) in type_map.items():
            if (param_type == c_type):
                param_type = cs_type
                break
        
        if param_type.startswith('const '):
            param_type = param_type.rsplit(' ', 1)[-1]
            
        params[i] = (param_type, param_name)
        
    return (ret_type, name, params, del_clazz)


def generate_interop(incl_file, ignore):
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

    functions = [convert((ret_type, name, params)) for (ret_type, name, params) in parsed if name not in ignore]

    functions_text = '\n\n'.join((get_function_text(f) for f in functions))

    interop_text = f'''using System.Runtime.InteropServices;

namespace Cantera.Interop;

static partial class LibCantera
{{
{functions_text}
}}'''

    with open(gen_file_dir + 'Interop.LibCantera.' + incl_file.name + '.g.cs', 'w') as f:
        f.write(interop_text)

    handles = [(name, del_clazz) for (_, name, _, del_clazz) in functions if del_clazz]

    if not handles:
        return

    handles_text = 'namespace Cantera.Interop;\n\n' + '\n\n'.join((get_base_handle_text(h) for h in handles))

    with open(gen_file_dir + 'Interop.Handles.' + incl_file.name + '.g.cs', 'w') as f:
        f.write(handles_text)


(_, include_dir, gen_file_dir) = sys.argv

with open('generate_interop.yaml', 'r') as config_file:
   config = ruamel.yaml.safe_load(config_file)

ignore = config['ignore']

for incl_file in os.scandir('Include/clib'):
    if incl_file.name not in ignore or ignore[incl_file.name]:
        generate_interop(incl_file, ignore.get(incl_file.name, []))

derived_handles = '\n\n'.join((get_derived_handle_text(d) for d in config['derived_handles'].items()))

derived_handles_text = f'''using System.Runtime.InteropServices;

namespace Cantera.Interop;

{derived_handles}'''

with open(gen_file_dir + 'Interop.Handles.g.cs', 'w') as f:
    f.write(derived_handles_text)