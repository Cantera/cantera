clear all

cantera_root = getenv('CANTERA_ROOT');
ctName = '/test/matlab_experimental/libcantera_shared.so';

% Unload Cantera and remove temporary library file
unloadlibrary('libcantera_shared');
delete([cantera_root, ctName]);

disp('Cantera has been unloaded');
