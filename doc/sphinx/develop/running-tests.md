# Running Tests and Debugging

## Using Cantera from the Build Directory

After compiling Cantera from source, it is convenient to use the Cantera Python module
directly without fully installing Cantera. Doing this requires some special
configuration to make sure your operating system and Python installation can find all
the relevant files.

### Temporary Configuration

To use the Python module from the build directory for the duration of the current
shell session, run the following commands from the root of the Cantera source directory:

::::{tab-set}
:::{tab-item} Linux
:sync: linux
```sh
export LD_LIBRARY_PATH=$(pwd)/build/lib
export PYTHONPATH=$(pwd)/build/python
```
:::

:::{tab-item} macOS
:sync: macos
```sh
export DYLD_LIBRARY_PATH=$(pwd)/build/lib
export PYTHONPATH=$(pwd)/build/python
```
:::

:::{tab-item} Windows (PowerShell)
:sync: windows
```pwsh
$Env:PYTHONPATH = (Get-Location).Path + '\build\python'
```
:::
::::

After this, you will be able to use the Cantera Python module from the terminal where
these commands were executed. For example:

```sh
# start an interactive shell and then import Cantera
ipython

# convert a mechanism
python -m cantera.ck2yaml --inp=...

# run an example
python samples/python/onedim/adiabatic_flame.py
```

### Permanent Configuration Using Conda

The following commands will make the Cantera Python module in the build directory
available within the current Conda environment.

::::{tab-set}
:::{tab-item} Linux
:sync: linux
```sh
conda env config vars set LD_LIBRARY_PATH=$(pwd)/build/lib
conda env config vars set PYTHONPATH=$(pwd)/build/python

```
:::

:::{tab-item} macOS
:sync: macos
```sh
conda env config vars set DYLD_LIBRARY_PATH=$(pwd)/build/lib
conda env config vars set PYTHONPATH=$(pwd)/build/python
```
:::

:::{tab-item} Windows (PowerShell)
:sync: windows
```pwsh
New-Item -Type Junction -Target ((Get-Location).Path + '\build\python\cantera') -Path ($Env:CONDA_PREFIX + '\lib\site-packages\cantera')
```
:::
::::

Setting up the paths this way can be useful if you're using an IDE such as VS Code that
is aware of Conda environments and may allow you to run tests within the IDE.

## Running Individual Python Tests and Test Subsets

To run the full Python test suite:

```sh
pytest -raP --verbose test/python
```

To run tests from a specific file:
```sh
pytest -raP --verbose test/python/test_composite.py
```

To run tests where the name matches a specific pattern:
```sh
pytest -raP --verbose test/python -k 'add_rxn'
```

To have the test runner activate a PDB shell after a failed assertion:
```sh
pytest -raP --verbose --pdb test/python```
```

```{important}
Remember to run `scons build` after any code change before running tests directly with
`pytest`.
```

## Running tests through GDB

If you are trying to debug a problem that results in a segmentation fault (crash)
during the test, you can try running it inside the debugger. The following notes may
help with running either the C++ or Python tests through GDB on Linux.

- Build the source code with debugging information enabled and optimizations disabled.
  It can also help to enable certain bounds checks available in the standard library:
  ```sh
  scons build optimize=n debug=y debug_flags='-g -D_GLIBCXX_ASSERTIONS'
  ```

### Running C++ (GoogleTest) Tests

First, run the target test through SCons to make sure the Cantera library and the test
executable are up to date. We'll use the `kinetics` test subset as an example, with the
goal of debugging the `PlogLowPressure` test.
```sh
scons test-kinetics
```

Next, we need to set some environment variables that specify where the test program
can find files needed to run the tests. Normally, these are handled by SCons:
```sh
export CANTERA_DATA=test/data:data
export PYTHONPATH=test/python:$PYTHONPATH
```

To start the debugger and enable only the desired test:
```sh
gdb --args build/test/kinetics/kinetics --gtest_filter="*PlogLowPressure"
```

Then at the `(gdb)` prompt, type `run` and press Enter.

After the debugger reaches the error, the most useful command is the `where` command,
which will print a stack trace. The `up` command will move up the stack each time it is
called and print the corresponding line of the source code.

