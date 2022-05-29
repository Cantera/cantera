from dataclasses import dataclass
import os
import typing

from .._types import * 
from .._types import _unpack


@dataclass
class _CsFunc:
    """ Represents a C# method """

    ret_type: str
    name: str
    params: list[Param]
    del_clazz: typing.Optional[str]

    
    def __iter__(self):
        return _unpack(self)


class SourceGenerator(SourceGeneratorBase):
    _prolog = '    [DllImport(LibFile)]\n    public static extern'

    _type_map = {
        'const char*': 'string',
        'size_t': 'nuint',
        'char*': 'byte*'
    }


    @staticmethod
    def _join_params(params):
        return ', '.join((p.p_type + ' ' + p.name for p in params))


    @classmethod
    def _get_function_text(cls, function):
        ref_type, name, params, _ = function
        is_unsafe = any((p.p_type.endswith('*') for p in params))
        if is_unsafe:
            return f'{cls._prolog} unsafe {ref_type} {name}({cls._join_params(params)});'
        else:
            return f'{cls._prolog} {ref_type} {name}({cls._join_params(params)});'


    @staticmethod
    def _get_base_handle_text(handle):
        name, del_clazz = handle

        handle = f'''class {del_clazz} : CanteraHandle
{{
    protected override bool ReleaseHandle() =>
        Convert.ToBoolean(LibCantera.{name}(Value));
}}'''

        return handle


    @staticmethod
    def _get_derived_handle_text(derived):
        derived, base = derived

        derived = f'''class {derived} : {base} {{ }}'''

        return derived


    @classmethod
    def _convert_func(cls, parsed: Func):
        ret_type, name, params = parsed
        clazz, method = name.split('_', 1)

        #copy the params list
        params = list(params)
        
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
                _, param_name = params[0]
                params[0] = handle_clazz, param_name
                
        for c_type, cs_type in cls._type_map.items():
            if ret_type == c_type:
                ret_type = cs_type
                break

        for i in range(0, len(params)):
            param_type, param_name = params[i]
            
            for c_type, cs_type in cls._type_map.items():
                if param_type == c_type:
                    param_type = cs_type
                    break
            
            if param_type.startswith('const '):
                param_type = param_type.rsplit(' ', 1)[-1]
                
            params[i] = Param(param_type, param_name)
            
        return _CsFunc(ret_type, name, params, del_clazz)


    def __init__(self, out_dir: str, config: dict):
        self._out_dir = out_dir
        self._config = config


    def generate_source(self, incl_file: os.DirEntry, funcs: list[Func]):
        cs_funcs = [self._convert_func(f) for f in funcs]

        functions_text = '\n\n'.join((self._get_function_text(f) for f in cs_funcs))

        interop_text = f'''using System.Runtime.InteropServices;

namespace Cantera.Interop;

static partial class LibCantera
{{
{functions_text}
}}'''

        with open(self._out_dir + 'Interop.LibCantera.' + incl_file.name + '.g.cs', 'w') as f:
            f.write(interop_text)

        handles = [(name, del_clazz) for _, name, _, del_clazz in cs_funcs if del_clazz]

        if not handles:
            return

        handles_text = 'namespace Cantera.Interop;\n\n' + '\n\n'.join((self._get_base_handle_text(h) for h in handles))

        with open(self._out_dir + 'Interop.Handles.' + incl_file.name + '.g.cs', 'w') as f:
            f.write(handles_text)


    def finalize(self):
        derived_handles = '\n\n'.join((self._get_derived_handle_text(d) for d in self._config['derived_handles'].items()))

        derived_handles_text = f'''using System.Runtime.InteropServices;

namespace Cantera.Interop;

{derived_handles}'''

        with open(self._out_dir + 'Interop.Handles.g.cs', 'w') as f:
            f.write(derived_handles_text)