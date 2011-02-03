The project setttings assume Cantera is installed in C:\CANTERA. 
If this is not the case, edit the project properties before building.

The project settings that differ from the defaults are:

1) use of multithreaded DLL system libraries
2) specification of the Cantera libraries as additional 
   input for the linker
3) specification of the Cantera include and lib directories.

These changes to the default project settings need to be made whenever 
you create a new project that you want to link to Cantera.