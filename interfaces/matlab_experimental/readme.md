# Experimental MATLAB Toolbox for Cantera
This experimental Matlab Toolbox for Cantera changes the Matlab interface to the modern
Matlab structure and syntaxes for OOP. It replaces the MEX interface with direct
function calling from Cantera CLib.

## Installation guide:

1. Install Matlab (any release newer than R2008a).
2. Compile Cantera from Source and install in your Conda environment, as directed in
   this link. https://cantera.org/stable/develop/compiling-install.html. The
   experimental Matlab Toolbox does not require a SCons option to install at this moment
   since it's stand-alone.
3. For first time users, launch Matlab, then navigate to `/path/to/cantera/source/code`
   (the folder containing `interfaces` and `samples`) using "Browse for Folder".
4. For Linux users: Matlab currently uses Intel MKL which uses 64-bit integer types
   that are incompatible with the standard 32-bit integers used by the default version
   of OpenBLAS that comes with Cantera. As such, the correct environment variables
   need to be set to launch Matlab with the correct BLAS/LAPACK libraries loaded: 
   `export LD_PRELOAD=/path/to/openblas/library:/path/to/lapack/library`
   then, launch Matlab in the terminal with:
   `matlab -softwareopengl`.
5. In the Matlab command window, run
   `addpath(genpath([pwd, '/interfaces/matlab_experimental']))` to add search path for
   the experimental toolbox.
6. In the Matlab command window, run
   `cd([pwd, '/interfaces/matlab_experimental/Utility'])` to navigate to the Utility
   folder.
7. Open the file named 'ctRoot.m', in the second line, edit `output=` to
   `output='/path/to/conda/environment'`, then save the file. This sets the search path
   for the `ctLoad` command to find the shared library file for Cantera.
8. In the Matlab command window, run `savepath` to save all search paths.
9. To switch back to the legacy Matlab toolbox, revert the search paths.

## Usage guide:

1. To start using the experimental toolbox, run `ctLoad` command.
2. Refer to examples in `/samples/matlab_experimental` in your
   Cantera source code directory.
3. To stop using the new Cantera interface, run the following commands:
   `ctCleanUp` `ctUnload`.
