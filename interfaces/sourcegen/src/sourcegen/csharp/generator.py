"""Generator for C# source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import sys
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from ..dataclasses import Func, Param, HeaderFile, ArgList
from ..generator import SourceGenerator

from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()
_LOADER = Environment(loader=BaseLoader)


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the CSharpSourceGenerator class."""

    c_type_crosswalk: dict[str, str]  #: C type crosswalks

    prop_type_crosswalk: dict[str, str]  #: Span type crosswalks

    class_crosswalk: dict[str, str]

    class_accessors: dict[str, str]

    handle_array_takers: dict[str, str]

    derived_handles: dict[str, str]

    wrapper_classes: dict[str, dict[str, str]]


@dataclass(frozen=True)
@with_unpack_iter
class CsFunc(Func):
    """Represents a C# interop method."""
    is_handle_release_func: bool = None
    handle_class_name: str | None = None

    def has_string_param(self) -> bool:
        """Identify any parameters that take strings."""
        return any(p.p_type == 'string' for p in self.arglist)

    def returns_handle(self) -> bool:
        """True if this function returns a handle."""
        return self.ret_type.endswith("Handle")

    def gets_string(self) -> bool:
        """True if this function is used to get a string."""
        return (len(self.arglist) >= 2
                and self.arglist[-1].p_type == "Span<byte>")

    def gets_double_array(self) -> bool:
        """True if this function is used to get an array of doubles."""
        return (len(self.arglist) >= 2
                and self.arglist[-1].p_type == "Span<double>")

    def sets_array(self) -> bool:
        """True if this function is used to set an array."""
        return (len(self.arglist) >= 2
                and self.arglist[-1].p_type.startswith("ReadOnlySpan"))


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
        getter_name = f"{clib_area}_{c_name}"
        getter = known_funcs.get(getter_name)
        if not getter:
            getter_name = f"{clib_area}_get{c_name.capitalize()}"
            getter = known_funcs.get(getter_name)
        if not getter:
            _LOGGER.critical(f"No getter found for {c_name!r}!")
            sys.exit(1)

        setter_name = f"{clib_area}_set{c_name.capitalize()}"
        setter = known_funcs.get(setter_name)
        if not setter:
            setter_name = None

        match len(getter.arglist):
            case 1:
                # This is a simple scalar property
                prop_type = getter.ret_type
            case 3:
                # This is a property that returns an array or string.
                # Note that _scaffold_interop generates wrappers that avoid needing
                # explicit size or buffer arguments. The property will use these
                # wrappers, not the function discovered here.
                prop_type = getter.arglist[-1].p_type
                prop_type = self._config.prop_type_crosswalk[prop_type]
            case _:
                _LOGGER.critical(f"Getter {getter_name!r} has an "
                    "unsupported signature!")
                sys.exit(1)

        if prop_type in ["int", "double", "string", "double[]"]:
            template = _LOADER.from_string(self._templates["csharp-property"])
            return template.render(
                prop_type=prop_type, cs_name=cs_name,
                getter=getter_name, setter=setter_name)

        _LOGGER.critical(f"Unable to scaffold properties of type {prop_type!r}!")
        sys.exit(1)

    def _get_wrapper_class_name(self, clib_area: str) -> str:
        return self._config.class_crosswalk[clib_area]

    def _get_handle_class_name(self, clib_area: str) -> str:
        return self._get_wrapper_class_name(clib_area) + "Handle"

    def _convert_func(self, parsed: Func) -> CsFunc:
        """Convert CLib signature to C# signature."""
        # TODO: The Func object contains information on CLib header and the underlying
        # C++ implementation. Some information (brief, wraps, returns, base, uses)
        # is currently preserved but not used.
        ret_type, name, params, brief, wraps, returns, base, uses = parsed

        clib_area, method = name.split("_", 1)

        # Shallow copy the params list
        # Some of the C# params will have the same syntax as the C params.
        # Others will be represented differently on the C# side, and we will
        # replace their entry in the list.
        # Therefore, copy the list so that we don’t accidentally modify
        # the params list which is attached to the C func.
        params = params[:]

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

        def crosswalk(par: str) -> str:
            """Crosswalk of C/C# types."""
            return self._config.c_type_crosswalk.get(par) or par

        ret_type = crosswalk(ret_type)

        handle_array_type = self._config.handle_array_takers.get(name)

        for i, param in enumerate(params):
            param_type = crosswalk(param.p_type)
            if handle_array_type and param_type.endswith("Span<int>"):
                # There is a slight inconsistency in CLib in that sometimes
                # the collection of handles will be marked as const and sometimes not.
                # However, is it never modified, so we can always use a ReadOnlySpan.
                # When this inconsistency is removed, the above can become
                # param_type == "ReadOnlySpan<int>"
                param_type = f"ReadOnlySpan<{handle_array_type}>"

            params[i] = Param(param_type, param.name, param.description,
                              param.direction, param.default, param.base)

        func = CsFunc(ret_type,
                      name,
                      ArgList(params),
                      brief,
                      wraps,
                      returns,
                      base,
                      uses,
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
        pinvoke_template = _LOADER.from_string(self._templates["csharp-interop-func"])
        function_list = [
            pinvoke_template.render(has_string_param=func.has_string_param(),
                                    declaration=func.declaration(),
                                    check_return=(not func.is_handle_release_func
                                                  and not func.returns_handle()),
                                    public=(not func.gets_string())
                                            and not func.gets_double_array()
                                            and not func.sets_array())
            for func in cs_funcs]

        # Add wrappers for functions that get strings.
        def transform_to_getstring_func(func: CsFunc) -> CsFunc:
            arglist = ArgList(func.arglist[:-2])
            return CsFunc('string', func.name, arglist, False, None)

        getstring_template = _LOADER.from_string(self._templates["csharp-getstring-func"])
        function_list += (
            getstring_template.render(declaration=(transform_to_getstring_func(func)
                                          .declaration()),
                                      invocation=func.invocation(),
                                      length_param_name=func.arglist[-2].name,
                                      span_param_name=func.arglist[-1].name)
            for func in cs_funcs if func.gets_string())

        # Add wrappers for functions that get or set arrays.
        def transform_to_span_func(func: CsFunc) -> CsFunc:
            arglist = ArgList([*func.arglist[:-2], func.arglist[-1]])
            return CsFunc('void', func.name, arglist, False, None)

        span_template = _LOADER.from_string(self._templates["csharp-span-func"])
        function_list += (
            span_template.render(declaration=(transform_to_span_func(func)
                                     .declaration()),
                                 invocation=func.invocation(),
                                 length_param_name=func.arglist[-2].name,
                                 span_param_name=func.arglist[-1].name)
            for func in cs_funcs if func.gets_double_array()
                                    or func.sets_array())

        # Add convenience overloads for functions that get arrays of doubles
        # to allocate the array and return it.
        def transform_to_getarray_func(func: CsFunc) -> CsFunc:
            arglist = ArgList([*func.arglist[:-2]])
            return CsFunc('double[]', func.name, arglist, False, None)

        def find_get_size_func(func: CsFunc) -> CsFunc:
            # Find the get-size function that wraps the C++ member
            # used by this function.
            return next(f for f in cs_funcs if f.wraps == func.uses[0])

        getarray_template = _LOADER.from_string(self._templates["csharp-getarray-func"])
        function_list += (
            getarray_template.render(declaration=(transform_to_getarray_func(func)
                                         .declaration()),
                                     get_size_invocation=(find_get_size_func(func)
                                         .invocation()),
                                     invocation=func.invocation(),
                                     length_param_name=func.arglist[-2].name,
                                     span_param_name=func.arglist[-1].name)
            for func in cs_funcs if func.gets_double_array() and func.uses)

        for f in cs_funcs:
            if f.gets_double_array() and not f.uses:
                _LOGGER.warning("Unable to generate simple get-array wrapper for "
                                f"{f.name!r} because the get-size function is unknown.")

        file_name = f"Interop.LibCantera.{header_file}.g.cs"
        self._write_file(
            file_name, "template_interop.g.cs.in", cs_functions=function_list)

    def _scaffold_handles(self, header_file: str, handles: dict[str, str]) -> None:
        template = _LOADER.from_string(self._templates["csharp-base-handle"])
        handle_list = [
            template.render(class_name=key, release_func_name=val)
            for key, val in handles.items()]

        file_name = f"Interop.Handles.{header_file}.g.cs"
        self._write_file(
            file_name, "template_handles.g.cs.in", cs_handles=handle_list)

    def _scaffold_derived_handles(self) -> None:
        template = _LOADER.from_string(self._templates["csharp-derived-handle"])
        handle_list = [
            template.render(derived_class_name=key, base_class_name=val)
            for key, val in self._config.derived_handles.items()]

        file_name = "Interop.Handles.g.cs"
        self._write_file(
            file_name, "template_handles.g.cs.in", cs_handles=handle_list)

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
        try:
            self._out_dir.mkdir(parents=True)
        except FileExistsError:
            # delete any existing files in the directory
            for f in self._out_dir.iterdir():
                f.unlink()

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
