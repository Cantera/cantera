(sec-dependencies)=
# Software used by Cantera

This section lists the versions of third-party software that are required to build and
use Cantera.

## Compilers

You must have one of the following C++ compilers installed on your system. A Fortran
compiler is required only if you plan to build the Fortran module.

- GNU compilers (C/C++/Fortran)

  - Known to work with versions 9.4 and 11.4. Expected to work with version >= 9.0.

- Clang/LLVM (C/C++)

  - Known to work with versions 10 and 12. Expected to work with version >= 5.0
  - Works with the versions included in Xcode 13.0 and 14.3.1.

- Intel compilers (C/C++/Fortran)

  - Known to work with the Intel OneAPI Compilers (version 2022.0.2).
  - Some earlier versions of the Intel compiler (including the 2017 version) are
    **NOT RECOMMENDED** because of a bug in the C compiler.

- Microsoft compilers (C/C++)

  - Known to work with Visual Studio 2019 (MSVC 14.2) and Visual Studio 2022
    (MSVC 14.3).

- MinGW (C/C++/Fortran)

  - <http://mingw-w64.org/doku.php> (64-bit)
  - <http://tdm-gcc.tdragon.net/> (64-bit)
  - Known to work with Mingw-w64 12.2.

## Other Required Software

- [SCons](https://scons.org/tag/releases.html):

  - Works with versions >= 4.5.0
  - On Windows, more recent SCons versions are required to support each new version of
    the MSVC compiler.

- [Python](https://python.org/downloads/):

  - Works with versions >= 3.12.

- [Boost](https://www.boost.org/releases/latest/)

  - Known to work with versions 1.71, 1.74, and 1.82; Expected to work with versions >=
    1.70
  - Only the "header-only" portions of Boost are required. Cantera does not currently
    depend on any of the compiled Boost libraries.

- [SUNDIALS](https://computing.llnl.gov/projects/sundials)

  - If SUNDIALS is not installed and you have checked out the Cantera source code using
    Git, SUNDIALS will be automatically downloaded and the necessary portions will be
    compiled and installed with Cantera.
  - Known to work with versions >= 6.4 and \<= 7.5.
  - To use SUNDIALS with Cantera on a Linux/Unix system, it must be compiled
    with the `-fPIC` flag. You can specify this flag when configuring SUNDIALS as
    `cmake -DCMAKE_C_FLAGS=-fPIC <other command-line options>`

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

  - If Eigen is not installed and you have checked out the Cantera source code using
    Git, Eigen will be automatically downloaded and installed with Cantera.
  - Known to work with version 3.4.0.

- [fmt](https://fmt.dev/latest/index.html)

  - If fmt (previously known as cppformat) is not installed and you have checked out the
    Cantera source code using Git, fmt will be automatically downloaded and the
    necessary portions will be compiled and installed with Cantera.
  - Known to work with versions 8.0 through 11.0.

- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

  - If yaml-cpp is not installed and you have checked out the Cantera source code using
    Git, it will be automatically downloaded and the necessary portions will be compiled
    and installed with Cantera.
  - Known to work with version 0.7.0. Version 0.6.0 or newer is required.

- [Doxygen](http://doxygen.nl/)

  - Required for building the C++ API Documentation
  - Required to generate the C and .NET interfaces.
  - Version 1.8 or newer is recommended.

- [Jinja2](https://jinja.palletsprojects.com/en/stable/)

  - Required for generated CLib and .NET.
  - Expected to work with versions >= 2.6.
  - Version 2.10 or newer is recommended.

## Optional Dependencies

- [NumPy](https://www.numpy.org/)

  - Required to build the Cantera Python module, and to run significant portions
    of the test suite.
  - Expected to work with versions >= 1.21.0.

- [Cython](https://cython.org/)

  - Required version >=3.0.8 to build the Python module.

- [pip](https://pip.pypa.io/en/stable/installing/) (Python)

  - Required to build the Cantera Python module.
  - Provides the `pip` command which can be used to install most of
    the other Python dependencies.

- [wheel](https://pypi.org/project/wheel/) (Python)

  - Required to build the Cantera Python module.

- [setuptools](https://pypi.org/project/setuptools/) (Python)

  - Required to build the Cantera Python module.

- [ruamel.yaml](https://pypi.org/project/ruamel.yaml/) (Python)

  - Required to convert input files from Chemkin, {term}`CTI`, and XML to the YAML
    format
  - Expected to work with versions >= 0.17.16

- [libhdf5](https://www.hdfgroup.org/solutions/hdf5/)

  - Required to read and write data files in the HDF5 format
  - Known to work with versions 1.12 and 1.14.

- [HighFive](https://highfive-devs.github.io/highfive/)

  - Required to read and write data files in the HDF5 format
  - If HighFive is not installed and you have checked out the Cantera source code
    using Git, HighFive will be automatically downloaded and the necessary portions will
    be compiled as part of the Cantera build process.
  - Version 2.5.0 or newer is required.

- [Google Test](https://github.com/google/googletest)

  - If Google Test is not installed and you have checked out the Cantera source code
    using Git, Google Test will be automatically downloaded and the necessary portions
    will be compiled as part of the Cantera build process.
  - Required to run significant portions of the test suite.
  - Known to work with version 1.11.0.

- [pytest](https://pytest.org)

  - Required to run the Python test suite.
  - Known to work with version 7.2.0

- MATLAB

  - Required to run the experimental MATLAB toolbox.
  - Known to work with R2025b. Expected to work with versions >= R2024a.

- [Windows Installer XML (WiX) toolset](http://wixtoolset.org/)

  - Required to build MSI installers on Windows.
  - Known to work with versions 3.5 and 3.8.

- Packages required for building Sphinx documentation

  - [Sphinx](https://www.sphinx-doc.org/en/stable/), version 7.0.0 or newer
  - [PyData Sphinx theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/),
    version 0.14.1 or newer
  - [MyST Parser](https://myst-parser.readthedocs.io/en/latest/), version 2.0.0 or newer
  - [MyST-NB](https://myst-nb.readthedocs.io/en/latest/), version 1.0.0 or newer
  - [Sphinx Gallery](https://sphinx-gallery.github.io/stable/index.html), version 0.15.0
    or newer
  - [Sphinx-copybutton](https://sphinx-copybutton.readthedocs.io/en/latest/), version
    0.5.2 or newer
  - [sphinx-design](https://sphinx-design.readthedocs.io/en/latest/), version 0.5.0 or
    newer
  - [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/),
    version 2.6.2 or newer
  - [sphinx-tags](https://github.com/Cantera/sphinx-tags)
    - Install forked version from Cantera organization on GitHub by running
      `pip install "git+https://github.com/Cantera/sphinx-tags.git@main"`
  - [Pygments](https://pygments.org/), version 2.17.2 or newer
  - [pyparsing](https://github.com/pyparsing/pyparsing/), version 3.1.2 or newer
  - [sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/), version 0.4.0
    or newer
  - [doxylink](https://pythonhosted.org/sphinxcontrib-doxylink/), version 1.12.3 or
    newer
  - [matlabdomain](https://pypi.org/project/sphinxcontrib-matlabdomain), version 0.21.5
    or newer
  - A working LaTeX installation (The conda `texlive-core` package will suffice)
  - Packages needed to run the examples for the example gallery:
    - [SciPy](https://scipy.org/)
    - [pandas](https://pandas.pydata.org/)
    - [CoolProp](http://www.coolprop.org/)
    - [Pint](https://pint.readthedocs.io/en/stable/)
    - [Matplotlib](https://matplotlib.org/)
    - [(python-)graphviz](https://graphviz.readthedocs.io/en/stable/)

- [Graphviz](https://www.graphviz.org/)

  - Required to build the dependency graph images in the C++ API Documentation
  - Known to work with version 2.40.1, expected to work with versions >=2.40.1

- [.NET](https://dotnet.microsoft.com/)

  - Required for the compilation of the experimental .NET interface
  - Known to work for the .NET 8.0 SDK
