# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import sys
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from ..dataclasses import Func, CFunc, Param, HeaderFile, ArgList
from ..generator import SourceGenerator

from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()
_LOADER = Environment(loader=BaseLoader)


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the CSharpSourceGenerator class."""

    ret_type_crosswalk: dict[str, str]  #: Return type crosswalks

    prop_type_crosswalk: dict[str, str]  #: Parameter type crosswalks

    class_crosswalk: dict[str, str]

    class_accessors: dict[str, str]

    derived_handles: dict[str, str]

    wrapper_classes: dict[str, dict[str, str]]


@dataclass(frozen=True)
@with_unpack_iter
class CsFunc(Func):
    """Represents a C# interop method"""
    # TODO: this should inherit from CFunc instead (or merged directly).
    # is_handle_release_func and handle_class_name may be inferred from existing
    # CFunc properties.

    is_handle_release_func: bool
    handle_class_name: str | None

    def has_string_param(self) -> bool:
        """Identify any parameters that take strings."""
        return any(p.p_type == 'string' for p in self.arglist)


class CSharpSourceGenerator(SourceGenerator):
    """The SourceGenerator for scaffolding C# files for the .NET interface"""

    def __init__(self, out_dir: str, config: dict, templates: dict) -> None:
        if not out_dir:
            _LOGGER.critical("Non-empty string identifying output path required.")
            sys.exit(1)
        self._out_dir = Path(out_dir)

        # use the typed config
        self._config = Config(**config)
        self._templates = templates

    def _get_property_text(self, clib_area: str, c_name: str, cs_name: str,
                           known_funcs: dict[str, CsFunc]) -> str:
        getter = known_funcs.get(clib_area + "3_" + c_name)

        if getter:
            # here we have found a simple scalar property
            prop_type = getter.ret_type
        else:
            # here we have found an array-like property (string, double[])
            getter = known_funcs[clib_area + "3_get" + c_name.capitalize()]
            # this assumes the last param in the function is a pointer type,
            # from which we determine the appropriate C# type
            prop_type = self._config.prop_type_crosswalk[getter.arglist[-1].p_type]

        setter = known_funcs.get(
            clib_area + "3_set" + c_name.capitalize(),
            CsFunc("", "", "", "", ""))

        if prop_type in ["int", "double"]:
            template = _LOADER.from_string(self._templates["csharp-property-int-double"])
            return template.render(
                prop_type=prop_type, cs_name=cs_name,
                getter=getter.name, setter=setter.name)

        if prop_type == "string":
            # for get-string type functions we need to look up the type of the second
            # (index 1) param for a cast because sometimes it's an int and other times
            # it's a nuint (size_t)
            template = _LOADER.from_string(self._templates["csharp-property-string"])
            return template.render(
                cs_name=cs_name, p_type=getter.arglist[1].p_type,
                getter=getter.name, setter=setter.name)

        _LOGGER.critical(f"Unable to scaffold properties of type {prop_type!r}!")
        sys.exit(1)

    def _get_wrapper_class_name(self, clib_area: str) -> str:
        return self._config.class_crosswalk[clib_area]

    def _get_handle_class_name(self, clib_area: str) -> str:
        return self._get_wrapper_class_name(clib_area) + "Handle"

    def _convert_func(self, parsed: CFunc) -> CsFunc:
        ret_type, name, _, _, _, _, _, _ = parsed
        clib_area, method = name.split("_", 1)
        clib_area = clib_area.rstrip("3")

        # Shallow copy the params list
        # Some of the C# params will have the same syntax as the C params.
        # Others will be represented differently on the C# side, and we will
        # replace their entry in the list.
        # Therefore, copy the list so that we don’t accidentally modify
        # the params list which is attached to the C func.
        params = parsed.arglist[:]

        release_func_handle_class_name = None

        if clib_area != "ct":
            handle_class_name = self._get_handle_class_name(clib_area)

            # It’s not a “global” function, therefore:
            #   * It wraps a constructor and returns a handle,
            #   * It wraps an instance method that returns a handle, or
            #   * It wraps an instance method and takes the handle as the first param.
            if method.startswith("del"):
                release_func_handle_class_name = handle_class_name
                params[0] = Param(handle_class_name, params[0].name)
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

        for i, param in enumerate(params):
            param_type = param.p_type
            param_name = param.name

            for c_type, cs_type in self._config.ret_type_crosswalk.items():
                if param_type == c_type:
                    param_type = cs_type
                    break

            # Most "setter" functions for arrays in CLib use a const double*,
            # but we also need to handle the cases for a plain double*
            if param_type == "double*" and method.startswith("set"):
                setter_double_arrays_count += 1
                if setter_double_arrays_count > 1:
                    # We assume a double* can reliably become a double[]. However, this
                    # logic is too simplistic if there is more than one array.
                    msg = f"Cannot scaffold {name!r} with multiple arrays of doubles!"
                    _LOGGER.critical(msg)
                    sys.exit(1)
                param_type = "double[]"

            params[i] = Param(param_type, param_name)

        func = CsFunc(ret_type,
                      name,
                      ArgList(params),
                      release_func_handle_class_name is not None,
                      release_func_handle_class_name)

        return func

    def _write_file(self, file_name: str, template_name: str, **kwargs) -> None:
        _LOGGER.info(f"  writing {file_name!r}")
        t_file = Path(__file__).parent / template_name
        template = _LOADER.from_string(t_file.read_text(encoding="utf-8"))
        contents = template.render(file_name=file_name, **kwargs)

        self._out_dir.joinpath(file_name).write_text(contents, encoding="utf-8")

    def _scaffold_interop(self, header_file: str, cs_funcs: list[CsFunc]) -> None:
        template = _LOADER.from_string(self._templates["csharp-interop-func"])
        function_list = [
            template.render(has_string_param=func.has_string_param(),
                            declaration=func.declaration())
            for func in cs_funcs]

        file_name = f"Interop.LibCantera.{header_file}.h.g.cs"
        self._write_file(
            file_name, "template_interop.h.g.cs.in", cs_functions=function_list)

    def _scaffold_handles(self, header_file: str, handles: dict[str, str]) -> None:
        template = _LOADER.from_string(self._templates["csharp-base-handle"])
        handle_list = [
            template.render(class_name=key, release_func_name=val)
            for key, val in handles.items()]

        file_name = f"Interop.Handles.{header_file}.h.g.cs"
        self._write_file(
            file_name, "template_handles.h.g.cs.in", cs_handles=handle_list)

    def _scaffold_derived_handles(self) -> None:
        template = _LOADER.from_string(self._templates["csharp-derived-handle"])
        handle_list = [
            template.render(derived_class_name=key, base_class_name=val)
            for key, val in self._config.derived_handles.items()]

        file_name = "Interop.Handles.g.cs"
        self._write_file(
            file_name, "template_handles.h.g.cs.in", cs_handles=handle_list)

    def _scaffold_wrapper_class(self, clib_area: str, props: dict[str, str],
                                known_funcs: dict[str, CsFunc]) -> None:
        property_list = [
            self._get_property_text(clib_area, c_name, cs_name, known_funcs)
            for c_name, cs_name in props.items()]

        wrapper_class_name = self._get_wrapper_class_name(clib_area)
        handle_class_name = self._get_handle_class_name(clib_area)
        file_name = wrapper_class_name + ".g.cs"
        self._write_file(
            file_name, "template_wrapper.g.cs.in",
            wrapper_class_name=wrapper_class_name, handle_class_name=handle_class_name,
            cs_properties=property_list)

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        self._out_dir.mkdir(parents=True, exist_ok=True)

        known_funcs: dict[str, list[CsFunc]] = {}

        for header_file in headers_files:
            cs_funcs = list(map(self._convert_func, header_file.funcs))
            known_funcs.update((f.name, f) for f in cs_funcs)

            file_name = header_file.path.name.replace("_auto.yaml", "")
            self._scaffold_interop(file_name, cs_funcs)

            handles = {func.handle_class_name: func.name
                for func in cs_funcs if func.is_handle_release_func}

            if not handles:
                continue

            self._scaffold_handles(file_name, handles)

        self._scaffold_derived_handles()

        for clib_area, props in self._config.wrapper_classes.items():
            self._scaffold_wrapper_class(clib_area, props, known_funcs)
