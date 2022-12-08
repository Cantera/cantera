This experimental Matlab Toolbox for Cantera changes the Matlab interface to the modern
Matlab structure and syntaxes for OOP. It replaces the MEX interface with direct
function calling from Canter Clib.

Installation guide:

1. Install Matlab (any release newer than R2008a).
2. Compile Cantera from Source and install in your Conda environment, as directed in
   this link. https://cantera.org/install/compiling-install.html. 2.5) The experimental
   Matlab Toolbox does not require a SCons option to install at this moment since it's
   stand-alone. It also does not require the current Matlab Toolbox to be installed.
3. For first time users, launch Matlab, then navigate to [/path/to/cantera/source/code]
   using "Browse for Folder". 3.5) Note for Ubuntu users, Matlab must be launched from
   the terminal using the following command:
   `LD_PRELOAD='/usr/lib/x86_64-linux-gnu/libstdc++.so.6' matlab`. This is because
   Matlab does not load the correct GLIBC++ library on start-up and will return an error
   when loading the Cantera toolbox.
4. In the Maltab command window, run
   `addpath(genpath([pwd, '/interfaces/matlab_experimental']))` to add search path for
   the experimental toolbox.
5. In the Maltab command window, run
   `addpath(genpath([pwd, '/samples/matlab_experimental']))` to add search path for the
   sample files.
6. In the Matlab command window, run
   `cd([pwd, '/interfaces/matlab_experimental/Utility'])` to navigate to the Utility
   folder.
7. Open the file named 'cantera_root.m', in the second line, edit `output=` to
   `output=[/path/to/conda/environment]`, then save the file. This sets the search path
   for the `LoadCantera` command to find the shared library file for Cantera.
8. After steps 3 to 7 are complete, make sure search paths to the current Matlab Toolbox
   and samples files are removed (if it's already installed). Having both the current
   and experimental version of the toolbox in the search path will lead to conflicts.
9. In the Matlab command window, run `savepath` to save all search paths.
10. To start using the experimental toolbox, run `LoadCantera` command.
11. To stop using the new Cantera interface, run the following commands: `clear all`
    `cleanup` `UnloadCantera`
12. To switch back to the current matlab toolbox, undo steps 3, 4, 5, 8, and 9. The
    command to remove search path in Matlab is `rmpath`.
13. A future updates will add automated installation.
