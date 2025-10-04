# Experimental MATLAB Toolbox for Cantera

The MATLAB Cantera interface has been fully modernized to improve usability,
performance, and maintainability. The legacy toolbox has been replaced with
a fully object-oriented structure, providing MATLAB classes that directly map to
Cantera objects and functions.

---

## Installation Guide

### Prerequisites

1. **MATLAB** (R2022a or later)
2. **MATLAB-compatible C++ compiler**. Refer to the [supported compilers](https://www.mathworks.com/support/requirements/supported-compilers.html) list.
3. **Conda package manager**

---

### Installation Steps

1. **Download the MATLAB Cantera toolbox**

   * Pull Cantera source code from [GitHub](https://github.com/Cantera/cantera)

2. **Obtain Cantera library and header files**

   * **Compile from source** and install in your Conda environment by following
   instructions: [Cantera Source Install](https://cantera.org/stable/develop/index.html)
   Note that the [`clib_legacy` configuration flag](https://cantera.org/3.1/develop/compiling/config-options.html) must be set to `true` for the current
   release.


3. **Add toolbox to MATLAB path**

   * **Windows/MacOS:** Launch MATLAB, navigate to the toolbox folder using *Browse for Folder*, and add it (including subfolders) to the MATLAB search path.
     Example:

     ```
     /path/to/cantera/interfaces/matlab_experimental
     ```
   * **Linux:** Configure MATLAB to use the system C++ library instead of MATLAB’s built-in version:

     ```bash
     export LD_PRELOAD=/path/to/system/C++/library
     ```

     Then launch MATLAB and add the toolbox path as above.

4. **Build the MATLAB Cantera interface**
   Run the following commands in MATLAB:

   ```matlab
   ctDir = '/path/to/cantera/source';
   ctIncludeDir = '/path/to/cantera/include';
   ctLibDir = '/path/to/cantera/library';
   ctBuildInterface(ctDir, ctIncludeDir, ctLibDir);
   ```

   * If you compiled Cantera from source, ctDir should be `path/to/cantera/include`;
   and ctLibDir should be `path/to/cantera/build/lib`.

5. **Verify build**
   After building, the compiled interface should appear under:

   ```
   /path/to/cantera/interfaces/matlab_experimental/cantera/ctMatlab
   ```

---

## Usage Guide

1. **Load the toolbox**

   ```matlab
   ctLoad
   ```

   The interface supports two modes:

   * `'outofprocess'` (default) — higher stability and compatibility; **required on Linux**.
   * `'inprocess'` — runs inside MATLAB with lower overhead.

2. **Run examples**
   Navigate to the samples folder:

   ```
   /samples/matlab_experimental
   ```

3. **Run the unit test suite**
   Navigate to the test folder:

   ```
   /test/matlab_experimental
   ```

   Execute:

   ```matlab
   runMatlabInterfaceTests.m
   ```

4. **Unload the toolbox**
   To stop using the interface, run:

   ```matlab
   ctCleanUp
   ctUnload
   ```

---
