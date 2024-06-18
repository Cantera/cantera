import sys
import os
import re
from Cython.Build import cythonize
from setuptools import setup, Extension
from setuptools.command.install import install
from setuptools.command.develop import develop
from pathlib import Path
import numpy
import shutil

PY_SRC = Path("cantera")
CT_SRC = Path("src")
EXT_SRC = Path("ext")
CT_INCLUDE = Path("include")
BOOST_INCLUDE = None
FORCE_CYTHON_COMPILE = False


class CanteraOptionsMixin:
    """Custom options for the install and develop commands.

    Modeled after https://stackoverflow.com/a/53833930
    """

    user_options = [
        ("force-cython-compile", None, "Force compilation of .pyx files via Cython"),
        ("boost-include", None, "Location of the Boost header files."),
    ]

    def initialize_options(self):
        super().initialize_options()
        self.force_cython_compile = False
        self.boost_include = None

    def finalize_options(self):
        if self.boost_include is not None:
            if not Path(self.boost_include).is_dir():
                raise TypeError(f"The path {self.boost_include!r} is not a directory.")
        super().finalize_options()

    def run(self):
        global BOOST_INCLUDE, FORCE_CYTHON_COMPILE
        BOOST_INCLUDE = self.boost_include
        FORCE_CYTHON_COMPILE = self.force_cython_compile
        super().run()


class InstallCommand(CanteraOptionsMixin, install):
    user_options = (
        getattr(install, "user_options", []) + CanteraOptionsMixin.user_options
    )


class DevelopCommand(CanteraOptionsMixin, develop):
    user_options = (
        getattr(develop, "user_options", []) + CanteraOptionsMixin.user_options
    )


pyx_sources = list(map(str, PY_SRC.glob("*.pyx")))
cythonize(
    pyx_sources,
    compiler_directives={"binding": True},
    force=FORCE_CYTHON_COMPILE,
)

ct_sources = list(map(str, CT_SRC.glob("**/*.cpp")))
cpp_sources = list(map(str, PY_SRC.glob("*.cpp")))
sundials_sources = list(map(str, EXT_SRC.glob("sundials/**/*.c")))
yaml_cpp_sources = list(map(str, EXT_SRC.glob("yaml-cpp/**/*.cpp")))
fmt_sources = list(map(str, EXT_SRC.glob("fmt/*.cc")))

include_dirs = [
    str(CT_INCLUDE),
    str(CT_INCLUDE / "cantera" / "ext"),
    str(CT_SRC),
    "cantera",
    numpy.get_include(),
]

if "BOOST_INCLUDE" in os.environ:
    include_dirs.append(os.environ["BOOST_INCLUDE"])
elif BOOST_INCLUDE is not None:
    include_dirs.append(BOOST_INCLUDE)


def configure_build():
    boost_version = ""
    boost_locs = (
        os.environ.get("BOOST_INCLUDE", None),
        BOOST_INCLUDE,
        "/usr/include",
        "/usr/local/include",
    )
    for boost_dir in boost_locs:
        if boost_dir is None:
            continue
        version_hpp = Path(boost_dir) / "boost" / "version.hpp"
        if not version_hpp.exists():
            continue
        boost_lib_version = re.search(
            r'^#define.*BOOST_LIB_VERSION.*"(\d+_\d+[_\d]*?)"$',
            version_hpp.read_text(),
            flags=re.MULTILINE,
        )
        if boost_lib_version is not None:
            boost_version = boost_lib_version.group(1)
            break

    if not boost_version:
        raise ValueError(
            "Could not find Boost headers. Please set an environment variable called "
            "BOOST_INCLUDE that contains the path to the Boost headers."
        )
    try:
        boost_minor_version = int(boost_version.split("_")[1])
    except ValueError:
        raise ValueError(
            f"Could not convert Boost minor version to integer: '{boost_version}'"
        ) from None
    if boost_minor_version < 61:
        raise ValueError("Cantera requires Boost version 1.61 or newer.")

    if sys.platform != "win32":
        extra_compile_flags = ["-std=c++17", "-g0"]
        sundials_configh = {
            "SUNDIALS_USE_GENERIC_MATH": "#define SUNDIALS_USE_GENERIC_MATH 1",
            "SUNDIALS_BLAS_LAPACK": "/* #undef SUNDIALS_BLAS_LAPACK */",
        }
        sundials_cflags = ["-w"]
        sundials_macros = []
    else:
        extra_compile_flags = ["/EHsc", "/std:c++17"]
        sundials_macros = [("_CRT_SECURE_NO_WARNINGS", None)]
        sundials_configh = {
            "SUNDIALS_USE_GENERIC_MATH": "/* #undef SUNDIALS_USE_GENERIC_MATH */",
            "SUNDIALS_BLAS_LAPACK": "/* #undef SUNDIALS_BLAS_LAPACK */",
        }
        sundials_cflags = []

    sun_config_h_in = Path("sundials_config.h.in").read_text()
    sun_config_h = Path("sundials_config.h")
    sun_config_h.write_text(sun_config_h_in.format_map(sundials_configh))
    shutil.copy2(sun_config_h, EXT_SRC / "sundials" / "sundials")
    shutil.copy2(sun_config_h, CT_INCLUDE / "cantera" / "ext" / "sundials")

    return extra_compile_flags, sundials_cflags, sundials_macros


if "bdist_wheel" in sys.argv:
    extra_compile_flags, sundials_cflags, sundials_macros = configure_build()
else:
    extra_compile_flags = []
    sundials_cflags = []
    sundials_macros = []


def lib_def(sources, cflags, include_dirs, macros):
    """Convenience factory to create the dictionary for a Setuptools library build."""
    return dict(
        sources=sources, cflags=cflags, include_dirs=include_dirs, macros=macros
    )


sundials_inc_dir = include_dirs + [str(EXT_SRC / "sundials" / "sundials")]
libraries = [
    (
        "sundials",
        lib_def(sundials_sources, sundials_cflags, sundials_inc_dir, sundials_macros),
    ),
    ("yaml-cpp", lib_def(yaml_cpp_sources, extra_compile_flags, include_dirs, [])),
    ("fmtlib", lib_def(fmt_sources, extra_compile_flags, include_dirs, [])),
]

extensions = [
    Extension(
        name="cantera._cantera",
        sources=cpp_sources + ct_sources,
        include_dirs=include_dirs,
        extra_compile_args=extra_compile_flags,
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    ),
]

setup(
    ext_modules=extensions,
    libraries=libraries,
    cmdclass={"install": InstallCommand, "develop": DevelopCommand},
)
