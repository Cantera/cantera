```{py:currentmodule} cantera
```

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
is aware of Conda environments and will help with running tests within the IDE.

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
pytest -raP --verbose --pdb test/python
```

```{important}
Remember to run `scons build` after any code change before running tests directly with
`pytest`.
```

(sec-pytest-vscode)=
## Running Python Tests Using VS Code

To allow running the Python unit tests through VS Code, add the following options
to `.vscode/settings.json`:
```json
      "python.testing.unittestEnabled": false,
      "python.testing.pytestEnabled": true,
      "python.testing.cwd": "test/python",
      "python.testing.unittestArgs": [
        "-v",
        "-s",
        "-p"
      ],
```

Once you have done this, you should be able to run and debug individual Python tests
directly from VS Code using the controls in the gutter to the left of the editor window.

(sec-gtest-vscode)=
## Debugging C++ Tests Using VS Code

To allow running the C++ unit tests through the [VS Code interactive
debugger](https://code.visualstudio.com/docs/editor/debugging), edit or create the file
`.vscode/launch.json` with contents based on the following template. Some substitutions
are needed for paths specific to your computer.

::::{tab-set}
:::{tab-item} Linux
:sync: linux
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "(gdb) Launch test/thermo",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/test/thermo/thermo",
            "cwd": "${workspaceFolder}/test/work",
            "environment": [
                {
                    "name": "LD_LIBRARY_PATH",
                    "value": "${workspaceFolder}/build/lib"
                },
                {
                    "name": "CANTERA_DATA",
                    "value": "${workspaceFolder}/data:${workspaceFolder}/test/data"
                },
                {
                    "name": "PYTHONPATH",
                    "value": "${workspaceFolder}/test/python:${workspaceFolder}/build/python:CONDAROOT/lib/python3.14/site-packages"
                }
            ],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        },
        {
            "name": "(gdb) Launch test/kinetics",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/test/kinetics/kinetics",
            // Optional configuration to run a specific test or subset of tests
            "args": [
                // "--gtest_filter=KineticsFromScratch.add_falloff_reaction"
            ],
            "cwd": "${workspaceFolder}/test/work",
            "environment": [
                {
                    "name": "LD_LIBRARY_PATH",
                    "value": "${workspaceFolder}/build/lib"
                },
                {
                    "name": "CANTERA_DATA",
                    "value": "${workspaceFolder}/data:${workspaceFolder}/test/data"
                },
                {
                    "name": "PYTHONPATH",
                    "value": "${workspaceFolder}/test/python:${workspaceFolder}/build/python:CONDAROOT/lib/python3.14/site-packages"
                }
            ],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        }
    ]
}
```
:::

:::{tab-item} macOS
:sync: macos
Install the [CodeLLDB](https://marketplace.visualstudio.com/items?itemName=vadimcn.vscode-lldb)
extension for VS Code.

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Thermo GTest",
            "type": "lldb",
            "request": "launch",
            "program": "${workspaceFolder}/build/test/thermo/thermo",
            "cwd": "${workspaceFolder}/test/work",
            "env": {
                "DYLD_LIBRARY_PATH": "${workspaceFolder}/build/lib",
                "CANTERA_DATA": "${workspaceFolder}/data:${workspaceFolder}/test/data",
                "PYTHONPATH": "${workspaceFolder}/test/python:${workspaceFolder}/build/python:CONDAROOT/lib/python3.14/site-packages"
            },
        },
        {
            "name": "Kinetics GTest",
            "type": "lldb",
            "request": "launch",
            "program": "${workspaceFolder}/build/test/kinetics/kinetics",
            // Optional configuration to run a specific test or subset of tests
            "args": [
                // "--gtest_filter=KineticsFromScratch.add_falloff_reaction"
            ],
            "cwd": "${workspaceFolder}/test/work",
            "env": {
                "DYLD_LIBRARY_PATH": "${workspaceFolder}/build/lib",
                "CANTERA_DATA": "${workspaceFolder}/data:${workspaceFolder}/test/data",
                "PYTHONPATH": "${workspaceFolder}/test/python:${workspaceFolder}/build/python:CONDAROOT/lib/python3.14/site-packages"
            },
        }
    ]
}
```
:::

:::{tab-item} Windows
:sync: windows
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "type": "cppvsdbg",
            "request": "launch",
            "name": "Thermo GTest",
            "program": "${workspaceFolder}/build/test/thermo/thermo.exe",
            "cwd": "${workspaceFolder}/test/work",
            "environment": [
                {
                    "name": "Path",
                    "value": "${env:Path};${workspaceFolder}/build/lib;CONDAROOT;CONDAROOT/Library/bin"
                },
                {
                    "name": "CANTERA_DATA",
                    "value": "${workspaceFolder}/data;${workspaceFolder}/test/data"
                },
                {
                    "name": "PYTHONPATH",
                    "value": "${workspaceFolder}/build/python;${workspaceFolder}/test/python"
                }
            ]
        },
        {
            "type": "cppvsdbg",
            "request": "launch",
            "name": "Kinetics GTest",
            "program": "${workspaceFolder}/build/test/kinetics/kinetics.exe",
            "cwd": "${workspaceFolder}/test/work",
            // Optional configuration to run a specific test or subset of tests
            "args": [
                // "--gtest_filter=KineticsFromScratch.add_falloff_reaction"
            ],
            "environment": [
                {
                    "name": "Path",
                    "value": "${env:Path};${workspaceFolder}/build/lib;CONDAROOT;CONDAROOT/Library/bin"
                },
                {
                    "name": "CANTERA_DATA",
                    "value": "${workspaceFolder}/data;${workspaceFolder}/test/data"
                },
                {
                    "name": "PYTHONPATH",
                    "value": "${workspaceFolder}/build/python;${workspaceFolder}/test/python"
                }
            ]
        }
    ]
}

```
:::
::::

In the above, replace `CONDAROOT` with your Conda environment's root directory. To
determine this path, run `conda env list`. For example, if you see output like
```
conda env list
# conda environments:
#
base                     C:\Users\speth\mambaforge
ctdev                 *  C:\Users\speth\mambaforge\envs\ctdev
```
then use `C:/Users/speth/mambaforge/envs/ctdev` as `CONDAROOT`, replacing forward
slashes with backslashes since the backslash is used as an escape character in JSON.

You will also need to replace any instances of `python3.14` in the above to refer to
the version of Python in your Conda environment.

You can create additional configurations for any of the test executables that are part
of the GTest test suite by changing the `program` name above.

:::{note}
Before running any of the GTest programs through VS Code, you will first need to
use SCons to compile the corresponding test program, for example by running
```sh
scons build-thermo
```

If you make any changes to the test program or the Cantera library code, re-run this
command before starting the debugger again to ensure that the changes take effect.
:::

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

First, build the target test through SCons to make sure the Cantera library and the test
executable are up to date. We'll use the `kinetics` test subset as an example, with the
goal of debugging the `PlogLowPressure` test.
```sh
scons build-kinetics
```

Next, you need to set some environment variables that specify where the test program can
find files needed to run the tests. Normally, setting these would be handled by SCons.
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

## Accessing Stack Traces

Before resorting to running in the debugger, it can sometimes be helpful to get the
stack trace associated with an error from within a normal workflow. To enable printing
a stacktrace after an unhandled exception or segmentation fault, you can call the C++
function {ct}`printStackTraceOnSegfault` or the Python function
{py:func}`print_stack_trace_on_segfault`. To add stack trace information to the error
message printed by any `CanteraError` exceptions, you can call the C++ method
{ct}`CanteraError::setStackTraceDepth` with a depth of 20 to print the top 20 frames of
the call stack, or similarly the Python method
{py:meth}`CanteraError.set_stack_trace_depth`.

## Diagnosing memory errors using Clang's AddressSanitizer

The Clang compiler includes the
[AddressSanitizer](https://clang.llvm.org/docs/AddressSanitizer.html) tool which
enables checking for a range of memory errors such as out-of-bounds array access and
using deleted/freed variables. To compile Cantera in debug mode with these checks
enabled, run SCons with the following options:

```sh
scons build CC=clang CXX=clang++ env_vars=all debug_flags="-g -fsanitize=address" debug_linker_flags="-lprofiler -fsanitize=address -shared-libasan"
```

On macOS, running the Python test suite also requires the following:
```
DYLD_INSERT_LIBRARIES=$(clang -print-file-name=libclang_rt.asan_osx_dynamic.dylib)
```

Or on Linux:
```
LD_PRELOAD=$(clang -print-file-name=libclang_rt.asan-x86_64.so)
```
