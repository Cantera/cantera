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

  - Known to work with Visual Studio 2017 (MSVC 14.1), Visual Studio 2019 (MSVC 14.2)
    and Visual Studio 2022 (MSVC 14.3).

- MinGW (C/C++/Fortran)

  - <http://mingw-w64.org/doku.php> (64-bit)
  - <http://tdm-gcc.tdragon.net/> (64-bit)
  - Known to work with Mingw-w64 12.2.

## Other Required Software

- SCons:

  - <http://scons.org/tag/releases.html>
  - Works with versions >= 3.0.0
  - On Windows, more recent SCons versions are required to support each new version of
    the MSVC compiler.

- Python:

  - <https://python.org/downloads/>
  - Works with versions >= 3.8.

- Boost

  - <https://www.boost.org/users/download/>
  - Known to work with versions 1.71, 1.74, and 1.82; Expected to work with versions >=
    1.70
  - Only the "header-only" portions of Boost are required. Cantera does not currently
    depend on any of the compiled Boost libraries.

- SUNDIALS

  - If SUNDIALS is not installed and you have checked out the Cantera source code using
    Git, SUNDIALS will be automatically downloaded and the necessary portions will be
    compiled and installed with Cantera.
  - <https://computing.llnl.gov/projects/sundials>
  - Known to work with versions >= 3.0 and \<= 7.0.
  - To use SUNDIALS with Cantera on a Linux/Unix system, it must be compiled
    with the `-fPIC` flag. You can specify this flag when configuring SUNDIALS as
    `cmake -DCMAKE_C_FLAGS=-fPIC <other command-line options>`

- Eigen

  - If Eigen is not installed and you have checked out the Cantera source code using
    Git, Eigen will be automatically downloaded and installed with Cantera.
  - <https://eigen.tuxfamily.org/index.php?title=Main_Page>
  - Known to work with version 3.4.0.

- fmt

  - If fmt (previously known as cppformat) is not installed and you have checked out the
    Cantera source code using Git, fmt will be automatically downloaded and the
    necessary portions will be compiled and installed with Cantera.
  - <https://fmt.dev/latest/index.html>
  - Known to work with version 9.1.0.

- yaml-cpp

  - If yaml-cpp is not installed and you have checked out the Cantera source code using
    Git, it will be automatically downloaded and the necessary portions will be compiled
    and installed with Cantera.
  - <https://github.com/jbeder/yaml-cpp>
  - Known to work with version 0.7.0. Version 0.6.0 or newer is required.

## Optional Dependencies

- [Numpy](https://www.numpy.org/)

  - Required to build the Cantera Python module, and to run significant portions
    of the test suite.
  - Expected to work with versions >= 1.12.0. 1.16.0 or newer is recommended.

- [Cython](https://cython.org/)

  - Required version >=0.29.31 to build the Python module.

- [Pip](https://pip.pypa.io/en/stable/installing/) (Python)

  - Required to build the Cantera Python module.
  - Provides the `pip` command which can be used to install most of
    the other Python dependencies.

- [wheel](https://pypi.org/project/wheel/) (Python)

  - Required to build the Cantera Python module.

- [setuptools](https://pypi.org/project/setuptools/) (Python)

  - Required to build the Cantera Python module.

- [Ruamel.yaml](https://pypi.org/project/ruamel.yaml/) (Python)

  - Required to convert input files from Chemkin, {term}`CTI`, and XML to the YAML
    format
  - Known to work with versions 0.15.42, 0.15.87, and 0.16.5
  - Expected to work with versions >= 0.15.0

- [libhdf5](https://www.hdfgroup.org/solutions/hdf5/)

  - Required to read and write data files in the HDF5 format
  - Known to work with versions 1.12 and 1.14.

- [HighFive](https://github.com/BlueBrain/HighFive)

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

- Matlab

  - Required to build the Cantera Matlab toolbox.
  - Known to work with 2022a. Expected to work with versions >= 2009a.

- [Windows Installer XML (WiX) toolset](http://wixtoolset.org/)

  - Required to build MSI installers on Windows.
  - Known to work with versions 3.5 and 3.8.

- Packages required for building Sphinx documentation

  - [Sphinx](http://www.sphinx-doc.org/en/stable/)
  - [PyData Sphinx theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/)
  - [Pygments](http://pygments.org/)
  - [pyparsing](https://sourceforge.net/projects/pyparsing/)
  - [sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/)
  - [doxylink](https://pythonhosted.org/sphinxcontrib-doxylink/)
  - [matlabdomain](https://pypi.org/project/sphinxcontrib-matlabdomain)

- [Doxygen](http://doxygen.nl/)

  - Required for building the C++ API Documentation
  - Version 1.8 or newer is recommended.

- [Graphviz](https://www.graphviz.org/)

  - Required to build the dependency graph images in the C++ API Documentation
  - Known to work with version 2.40.1, expected to work with versions >=2.40.1
