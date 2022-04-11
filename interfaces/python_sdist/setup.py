import sys
import os
import re
from setuptools import setup, Extension
from setuptools.command.install import install
from setuptools.command.develop import develop
from pathlib import Path
import numpy
import shutil

HERE = Path(__file__).parent
CT_SRC = HERE / "src"
EXT_SRC = HERE / "ext"
CT_INCLUDE = HERE / "include"
BOOST_INCLUDE = None
FORCE_CYTHON_COMPILE = False

CYTHON_BUILT_FILES = [HERE / "cantera" / f"_cantera.{ext}" for ext in ("cpp", "h")]


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
    user_options = (getattr(install, "user_options", [])
                    + CanteraOptionsMixin.user_options)


class DevelopCommand(CanteraOptionsMixin, develop):
    user_options = (getattr(develop, "user_options", [])
                    + CanteraOptionsMixin.user_options)


if (
    not all(p.exists() for p in CYTHON_BUILT_FILES)
    or "sdist" in sys.argv
    or FORCE_CYTHON_COMPILE
    or os.environ.get("FORCE_CYTHON_COMPILE", False)
):
    from Cython.Build import cythonize
    CYTHON_EXT = ".pyx"
    for p in CYTHON_BUILT_FILES:
        if p.exists():
            p.unlink()
else:
    CYTHON_EXT = ".cpp"

    def cythonize(extensions):
        """Define a no-op for when we're not using Cython."""
        return extensions

source_files = ["cantera/_cantera" + CYTHON_EXT]
source_files += list(map(str, CT_SRC.glob("**/*.cpp")))
sundials_sources = list(map(str, EXT_SRC.glob("sundials/**/*.c")))
yaml_cpp_sources = list(map(str, EXT_SRC.glob("yaml-cpp/**/*.cpp")))
fmt_sources = list(map(str, EXT_SRC.glob("fmt/*.cc")))
libexecstream_sources = [str(EXT_SRC / "libexecstream" / "exec-stream.cpp")]

include_dirs = [
    str(CT_INCLUDE),
    str(CT_INCLUDE / "cantera" / "ext"),
    str(CT_SRC),
    numpy.get_include()
]

if "BOOST_INCLUDE" in os.environ:
    include_dirs.append(os.environ["BOOST_INCLUDE"])
elif BOOST_INCLUDE is not None:
    include_dirs.append(BOOST_INCLUDE)


def configure_build():
    boost_version = ""
    boost_locs = (os.environ.get("BOOST_INCLUDE", None), BOOST_INCLUDE, "/usr/include",
                  "/usr/local/include")
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

    config_h = {}

    def create_config(key, value):
        config_h[key] = f"#define {key} {value}"

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
    if boost_minor_version < 56:
        create_config("CT_USE_DEMANGLE", 0)
    else:
        create_config("CT_USE_DEMANGLE", 1)

    if sys.platform == "darwin":
        create_config("DARWIN", 1)
    else:
        create_config("DARWIN", 0)

    # I have no idea if this is the value for Solaris, and I don't have a Solaris
    # machine handy to test. YOLO!
    if sys.platform.startswith("sunos"):
        create_config("SOLARIS", 1)
    else:
        create_config("SOLARIS", 0)

    if sys.platform != "win32":
        extra_compile_flags = ["-std=c++11", "-g0"]
        sundials_configh = {
            "SUNDIALS_USE_GENERIC_MATH": "#define SUNDIALS_USE_GENERIC_MATH 1",
            "SUNDIALS_BLAS_LAPACK": "/* #undef SUNDIALS_BLAS_LAPACK */"
        }
        sundials_cflags = ["-w"]
        sundials_macros = []
    else:
        extra_compile_flags = []
        sundials_macros = [("_CRT_SECURE_NO_WARNINGS", None)]
        sundials_configh = {
            "SUNDIALS_USE_GENERIC_MATH": "/* #undef SUNDIALS_USE_GENERIC_MATH */",
            "SUNDIALS_BLAS_LAPACK": "/* #undef SUNDIALS_BLAS_LAPACK */"
        }
        sundials_cflags = []

    sun_config_h_in = (HERE / "sundials_config.h.in").read_text()
    sun_config_h = HERE / "sundials_config.h"
    sun_config_h.write_text(sun_config_h_in.format_map(sundials_configh))
    shutil.copy2(sun_config_h, EXT_SRC / "sundials" / "sundials")
    shutil.copy2(sun_config_h, CT_INCLUDE / "cantera" / "ext" / "sundials")

    config_h_in = (HERE / "config.h.in").read_text()
    ct_config_h = HERE / "include" / "cantera" / "base" / "config.h"
    ct_config_h.write_text(config_h_in.format_map(config_h))

    return extra_compile_flags, sundials_cflags, sundials_macros


if "bdist_wheel" in sys.argv:
    extra_compile_flags, sundials_cflags, sundials_macros = configure_build()
else:
    extra_compile_flags = []
    sundials_cflags = []
    sundials_macros = []

extensions = cythonize([
    Extension(
        "cantera._cantera",
        source_files,
        include_dirs=include_dirs,
        extra_compile_args=extra_compile_flags,
        language="c++",
    ),
])


def lib_def(sources, cflags, include_dirs, macros):
    """Convenience factory to create the dictionary for a Setuptools library build."""
    return dict(sources=sources, cflags=cflags, include_dirs=include_dirs,
                macros=macros)


sundials_inc_dir = include_dirs + [str(EXT_SRC / "sundials" / "sundials")]
libraries = [
    ("sundials", lib_def(sundials_sources, sundials_cflags, sundials_inc_dir,
                         sundials_macros)),
    ("yaml-cpp", lib_def(yaml_cpp_sources, extra_compile_flags, include_dirs, [])),
    ("fmtlib", lib_def(fmt_sources, extra_compile_flags, include_dirs, [])),
    ("libexecstream", {"sources": libexecstream_sources, "include_dirs": include_dirs}),
]

setup(
    ext_modules=extensions,
    libraries=libraries,
    cmdclass={"install": InstallCommand, "develop": DevelopCommand},
)
