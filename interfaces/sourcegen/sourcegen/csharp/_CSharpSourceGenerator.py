# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from itertools import starmap
from pathlib import Path
from typing import List, Dict
import re

from ._dataclasses import CsFunc
from ._Config import Config
from .._helpers import normalize_indent
from .._dataclasses import Func, Param, HeaderFile
from .._SourceGenerator import SourceGenerator


class CSharpSourceGenerator(SourceGenerator):
    """The SourceGenerator for scaffolding C# files for the .NET interface"""

    @staticmethod
    def _join_params(params: List[Param]) -> str:
        return ", ".join(p.p_type + " " + p.name for p in params)

    def _get_interop_func_text(self, func: CsFunc) -> str:
        ret_type, name, params, _, _ = func
        requires_unsafe_keyword = any(p.p_type.endswith("*") for p in params)
        params_text = self._join_params(params)

        if requires_unsafe_keyword:
            return f"{self._config.func_prolog} unsafe {ret_type} {name}({params_text});"
        else:
            return f"{self._config.func_prolog} {ret_type} {name}({params_text});"

    @staticmethod
    def _get_base_handle_text(class_name: str, release_func_name: str) -> str:
        handle = normalize_indent(f"""
            class {class_name} : CanteraHandle
            {{
                protected override bool ReleaseHandle() =>
                    LibCantera.{release_func_name}(Value) == InteropConsts.Success;
            }}
        """)

        return handle

    @staticmethod
    def _get_derived_handle_text(derived_class_name: str, base_class_name: str) -> str:
        derived_text = f"""class {derived_class_name} : {base_class_name} {{ }}"""

        return derived_text


    def _get_property_text(self, clib_area: str, c_name: str, cs_name: str,
                           known_funcs: Dict[str, CsFunc]) -> str:
        getter = known_funcs.get(clib_area + "_" + c_name)

        if getter:
            # here we have found a simple scalar property
            prop_type = getter.ret_type
        else:
            # here we have found an array-like property (string, double[])
            getter = known_funcs[clib_area + "_get" + c_name.capitalize()]
            # this assumes the last param in the function is a pointer type,
            # from which we determine the appropriate C# type
            prop_type = self._config.prop_type_crosswalk[getter.params[-1].p_type]

        setter = known_funcs.get(clib_area + "_set" + c_name.capitalize())

        if prop_type in ["int", "double"]:
            text = f"""
                public {prop_type} {cs_name}
                {{
                    get => InteropUtil.CheckReturn(
                        LibCantera.{getter.name}(_handle));"""

            if setter:
                text += f"""
                    set => InteropUtil.CheckReturn(
                        LibCantera.{setter.name}(_handle, value));"""

            text += """
                }
            """
        elif prop_type == "string":
            p_type = getter.params[1].p_type

            # for get-string type functions we need to look up the type of the second
            # (index 1) param for a cast because sometimes it"s an int and other times
            # its a nuint (size_t)
            text = f"""
                public unsafe string {cs_name}
                {{
                    get => InteropUtil.GetString(40, (length, buffer) =>
                        LibCantera.{getter.name}(_handle, ({p_type}) length, buffer));
            """

            if setter:
                text += f"""
                    set => InteropUtil.CheckReturn(
                        LibCantera.{setter.name}(_handle, value));"""

            text += """
                }
            """
        else:
            raise ValueError(f"Unable to scaffold properties of type {prop_type}!")

        return normalize_indent(text)

    def __init__(self, out_dir: Path, config: dict):
        self._out_dir = out_dir

        # use the typed config
        self._config = Config.from_parsed(config)

    def _get_wrapper_class_name(self, clib_area: str) -> str:
        return self._config.class_crosswalk[clib_area]

    def _get_handle_class_name(self, clib_area: str) -> str:
        return self._get_wrapper_class_name(clib_area) + "Handle"

    def _convert_func(self, parsed: Func) -> CsFunc:
        ret_type, name, _ = parsed
        clib_area, method = name.split("_", 1)

        # Shallow copy the params list
        # Some of the C# params will have the same syntax as the C params.
        # Others will be represented differently on the C# side, and we will
        # replace their entry in the list.
        # Therefore, copy the list so that we don’t accidentally modify
        # the params list which is attached to the C func.
        params = parsed.params[:]

        release_func_handle_class_name = None

        if clib_area != "ct":
            handle_class_name = self._get_handle_class_name(clib_area)

            # It’s not a “global” function, therefore:
            #   * It wraps a constructor and returns a handle,
            #   * It wraps an instance method that returns a handle, or
            #   * It wraps an instance method and takes the handle as the first param.
            if method.startswith("del"):
                release_func_handle_class_name = handle_class_name
            elif method.startswith("new"):
                ret_type = handle_class_name
            elif name in self._config.class_accessors:
                ret_type = self._config.class_accessors[name]
                params[0] = Param(handle_class_name, params[0].name)
            elif params:
                params[0] = Param(handle_class_name, params[0].name)

        for c_type, cs_type in self._config.ret_type_crosswalk.items():
            if ret_type == c_type:
                ret_type = cs_type
                break

        setter_double_arrays_count = 0

        for i, (param_type, param_name) in enumerate(params):

            for c_type, cs_type in self._config.ret_type_crosswalk.items():
                if param_type == c_type:
                    param_type = cs_type
                    break

            # Most "setter" functions for arrays in CLib use a const double*,
            # but we also need to handle the cases for a plain double*
            if param_type == "double*" and method.startswith("set"):
                setter_double_arrays_count += 1
                if setter_double_arrays_count > 1:
                    # We assume a double* can reliably become a double[].
                    # However, this logic is too simplistic if there is
                    # more than one array.
                    raise ValueError(f"Cannot scaffold {name} with "
                        + "more than one array of doubles!")

                if clib_area == "thermo" and re.match("^set_[A-Z]{2}$", method):
                    # Special case for the functions that set thermo pairs
                    # This allows the C# side to pass a pointer to the stack
                    # Rather than allocating an array on the heap (which requires GC)
                    param_type = "(double, double)*"
                else:
                    param_type = "double[]"

            params[i] = Param(param_type, param_name)

        func = CsFunc(ret_type,
                      name,
                      params,
                      release_func_handle_class_name is not None,
                      release_func_handle_class_name)

        return func

    def _write_file(self, filename: str, contents: str):
        print("  writing " + filename)

        self._out_dir.joinpath(filename).write_text(contents)

    def _scaffold_interop(self, header_file_path: Path, cs_funcs: List[CsFunc]):
        functions_text = "\n\n".join(map(self._get_interop_func_text, cs_funcs))

        interop_text = normalize_indent(f"""
            {normalize_indent(self._config.preamble)}

            using System.Runtime.InteropServices;

            namespace Cantera.Interop;

            static partial class LibCantera
            {{
                {normalize_indent(functions_text)}
            }}
        """)

        self._write_file("Interop.LibCantera." + header_file_path.name + ".g.cs",
            interop_text)

    def _scaffold_handles(self, header_file_path: Path, handles: Dict[str, str]):
        handles_text = "\n\n".join(starmap(self._get_base_handle_text, handles.items()))

        handles_text = normalize_indent(f"""
            {normalize_indent(self._config.preamble)}

            namespace Cantera.Interop;

            {normalize_indent(handles_text)}
        """)

        self._write_file("Interop.Handles." + header_file_path.name + ".g.cs",
            handles_text)

    def _scaffold_derived_handles(self):
        derived_handles = "\n\n".join(starmap(self._get_derived_handle_text,
            self._config.derived_handles.items()))

        derived_handles_text = normalize_indent(f"""
            {normalize_indent(self._config.preamble)}

            namespace Cantera.Interop;

            {derived_handles}
        """)

        self._write_file("Interop.Handles.g.cs", derived_handles_text)

    def _scaffold_wrapper_class(self, clib_area: str, props: Dict[str, str],
                                known_funcs: Dict[str, CsFunc]):
        wrapper_class_name = self._get_wrapper_class_name(clib_area)
        handle_class_name = self._get_handle_class_name(clib_area)

        properties_text = "\n\n".join(
            self._get_property_text(clib_area, c_name, cs_name, known_funcs)
                for (c_name, cs_name) in props.items())

        wrapper_class_text = normalize_indent(f"""
            {normalize_indent(self._config.preamble)}

            using Cantera.Interop;

            namespace Cantera;

            public partial class {wrapper_class_name} : IDisposable
            {{
                readonly {handle_class_name} _handle;

                #pragma warning disable CS1591

                {normalize_indent(properties_text)}

                #pragma warning restore CS1591

                /// <summary>
                /// Frees the underlying resources used by the
                /// native Cantera library for this instance.
                /// </summary>
                public void Dispose() =>
                    _handle.Dispose();
            }}
        """)

        self._write_file(wrapper_class_name + ".g.cs", wrapper_class_text)

    def generate_source(self, headers_files: List[HeaderFile]):
        self._out_dir.mkdir(parents=True, exist_ok=True)

        known_funcs: Dict[str, List[CsFunc]] = {}

        for header_file in headers_files:
            cs_funcs = list(map(self._convert_func, header_file.funcs))
            known_funcs.update((f.name, f) for f in cs_funcs)

            self._scaffold_interop(header_file.path, cs_funcs)

            handles = {func.handle_class_name: func.name
                for func in cs_funcs if func.is_handle_release_func}

            if not handles:
                continue

            self._scaffold_handles(header_file.path, handles)

        self._scaffold_derived_handles()

        for (clib_area, props) in self._config.wrapper_classes.items():
            self._scaffold_wrapper_class(clib_area, props, known_funcs)
