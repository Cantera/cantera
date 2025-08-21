# Experimental MATLAB Toolbox for Cantera
This experimental Matlab Toolbox for Cantera changes the Matlab interface to the modern
Matlab structure and syntaxes for OOP. It replaces the MEX interface with direct
function calling from Cantera CLib.

## Installation guide:

1. Install Matlab (any release newer than R2008a), then choose one of the following
   methods to obtain the Cantera CLib and header files for the toolbox:
2. Compile Cantera from Source and install in your Conda environment, as directed in
   this link: https://cantera.org/stable/develop/index.html. OR
3. Install the Cantera development interface, as instructed in this link:
   https://cantera.org/stable/install/conda.html#development-c-fortran-90-interface.
4. For first time users, launch Matlab, then navigate to
   `/path/to/cantera/matlab/toolbox` using "Browse for Folder".
5. For Linux users: Matlab currently uses Intel MKL which uses 64-bit integer types
   that are incompatible with the standard 32-bit integers used by the default version
   of OpenBLAS that comes with Cantera. As such, the correct environment variables
   need to be set to launch Matlab with the correct BLAS/LAPACK libraries loaded:
   `export LD_PRELOAD=/path/to/openblas/library:/path/to/lapack/library`
   then, launch Matlab in the terminal with:
   `matlab -softwareopengl`.
6. Run the built-in Utility function `ctPaths('Set', envPath)` to set search paths
   for the Toolbox:
   - `envPath`: this is the location of the Conda environment where the CLib and header
   files are located.

## Usage guide:

1. To start using the experimental toolbox, run `ctLoad` command.
2. Refer to examples in `/samples/matlab_experimental` in your
   Cantera source code directory.
3. To run the unit test suite, navigate to `/test/matlab_experimental`,
   then run the script `runMatlabInterfaceTests.m`.
4. To stop using the new Cantera interface, run the following commands:
   `ctCleanUp`
   `ctUnload`.
5. To remove the Cantera paths from MATLAB preferences, run `ctPaths('Unset')`.
