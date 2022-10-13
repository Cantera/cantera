This new Matlab Toolbox for Cantera changed the Matlab interface to the modern Matlab structure and syntaxes for OOP. It replaced the MEX interface with direct function calling from Canter Clib.

How to install the new Toolbox (for now):

1) Install Matlab (any version newer than R2008a). The toolbox have been tested for R2020a, R2020b, R2021a, and R2021b.
2) Compile Cantera from Source, as directed in this link.
https://cantera.org/install/compiling-install.html
3) Here is a list of cantera.conf options I used for development:
python_package = 'full'
matlab_toolbox = 'y'
4) After building and installing Cantera from source. Move [path/to/cantera/source/code/interfaces/Matlab_Toolbox_Revamp] to [/path/to/cantera/installation/matlab]
5) Launch Matlab, then navigate to [/path/to/cantera/installation/matlab/Matlab_Toolbox_Revamp] using "Browse for Folder".
6) In the Matlab command window, type SetCanteraPath('path/to/cantera/installation'). This should generate a cantera_root.mat file under the Matlab_Toolbox_Revamp/Utility folder, which stores the path as a Matlab variable.
7) In the command window, type LoadCantera to start using Cantera.
8) The Example folder contains test examples.
9) To stop using Cantera, type the following commands:
>> clear all
>> close all
>> cleanup
>> UnloadCantera
