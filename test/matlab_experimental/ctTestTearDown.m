clear all

rootDir = fullfile(pwd);
ctName = '/test/matlab_experimental/libcantera_shared.so';

% Unload Cantera and remove temporary library file
unloadlibrary('libcantera_shared');
delete([rootDir, ctName]);

disp('Cantera has been unloaded');
