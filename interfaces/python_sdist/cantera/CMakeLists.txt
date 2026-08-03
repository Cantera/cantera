file(GLOB_RECURSE CY_SOURCES "*.pyx")
set(CY_OUTPUTS)
foreach(CY_SOURCE IN LISTS CY_SOURCES)
  cython_transpile(
    ${CY_SOURCE}
    LANGUAGE CXX
    CYTHON_ARGS --directive binding=True -X language_level=3
    OUTPUT_VARIABLE CPP_OUTPUT
    PUBLIC_HEADER_VARIABLE DELEGATOR_H
  )
  list(APPEND CY_OUTPUTS ${CPP_OUTPUT})
  if(${CY_SOURCE} MATCHES "_delegate_callbacks.pyx$")
    add_custom_target(DELEGATOR_H_EXT ALL DEPENDS ${DELEGATOR_H})
  endif()
endforeach()

# The aggregator module _cantera uses pure-Python (.py) Cython syntax but is always
# compiled INTO the merged _cantera extension (it hosts the metapath finder and the
# PyInit functions for the .pyx modules merged here), so cythonize it alongside the .pyx
# sources rather than as a standalone pure-Python module.
# set(CANTERA_AGG_CPP "${CMAKE_CURRENT_BINARY_DIR}/_cantera.cpp")
cython_transpile(
  _cantera.py
  LANGUAGE CXX
  CYTHON_ARGS --directive binding=True -X language_level=3
  OUTPUT_VARIABLE CANTERA_AGG_CPP
)
list(APPEND CY_OUTPUTS ${CANTERA_AGG_CPP})

# Cythonize each Cython pure-Python-syntax modules. How they are compiled depends on the
# platform (see below).
set(PURE_CY_MODULES _utils jacobians constants units func1 reactionpath
                    yamlwriter speciesthermo mixture transport solutionbase
                    kinetics reaction thermo reactor _onedim delegator)
set(PURE_CY_OUTPUTS)
foreach(MODULE IN LISTS PURE_CY_MODULES)
  cython_transpile(
    "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.py"
    LANGUAGE CXX
    CYTHON_ARGS --directive binding=True -X language_level=3
    OUTPUT_VARIABLE CPP_OUTPUT
  )
  list(APPEND PURE_CY_OUTPUTS ${CPP_OUTPUT})
endforeach()

if(CANTERA_PYODIDE)
  # On Pyodide the standalone-extension + shared-libcantera layout isn't usable,
  # so compile the pure-Python modules INTO the merged _cantera extension.
  # cantera.<module> is then provided by PyInit_<module> within _cantera and
  # resolved by the CythonPackageMetaPathFinder in _cantera.py; the raw .py is excluded
  # from the Pyodide wheel (see pyproject.toml.in) so it doesn't shadow the merged
  # module on import.
  list(APPEND CY_OUTPUTS ${PURE_CY_OUTPUTS})
endif()

# Give an extension an RPATH so it finds the shared libcantera beside itself in the
# package directory. (Not used on Pyodide, where libcantera is static and merged in.)
function(set_extension_rpath target)
  if(APPLE)
    set_target_properties(${target} PROPERTIES INSTALL_RPATH "@loader_path")
  elseif(NOT WIN32)
    set_target_properties(${target} PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()
endfunction()

# The extension manager is compiled here rather than with the rest of the
# Cantera source code because it needs _delegate_callbacks.h, built by the Cython step
# above.
add_library(ext_manager ../src/extensions/PythonExtensionManager.cpp)
target_include_directories(ext_manager PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}"
  "../include"
  "${Python_INCLUDE_DIRS}")
target_link_libraries(ext_manager PRIVATE fmt::fmt Eigen3::Eigen)
add_dependencies(ext_manager DELEGATOR_H_EXT)

python_add_library(_cantera MODULE ${CY_OUTPUTS} WITH_SOABI)
target_include_directories(_cantera PRIVATE "${Python_NumPy_INCLUDE_DIRS}" "../include")
target_compile_definitions(_cantera PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(_cantera PRIVATE ext_manager cantera_lib)
if(NOT CANTERA_PYODIDE)
  set_extension_rpath(_cantera)

  # On native platforms, build each pure-Python Cython module as its own extension linked
  # to the shared libcantera. The compiled module is installed into the package; its .py
  # source ships alongside it (auto-included by scikit-build) and is shadowed at import
  # time by the extension. (On Pyodide these modules are merged into _cantera instead --
  # see the cythonize step above.)
  foreach(MODULE IN LISTS PURE_CY_OUTPUTS)
    get_filename_component(MODULE_TARGET "${MODULE}" NAME_WE)
    python_add_library(${MODULE_TARGET} MODULE ${MODULE} WITH_SOABI)
    target_include_directories(${MODULE_TARGET} PRIVATE "${Python_NumPy_INCLUDE_DIRS}" "../include")
    target_compile_definitions(${MODULE_TARGET} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_link_libraries(${MODULE_TARGET} PRIVATE cantera_lib)
    set_extension_rpath(${MODULE_TARGET})
    install(TARGETS ${MODULE_TARGET} DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
  endforeach()
endif()
