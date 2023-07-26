clear all

cantera_root = getenv('CANTERA_ROOT');

if ispc
    ctName = '/test/matlab_experimental/cantera_shared.dll';
elseif ismac
    ctname = '/test/matlab_experimental/libcantera_shared.dylib';
elseif isunix
    ctname = '/test/matlab_experimental/libcantera_shared.so';
end

% Unload Cantera and remove temporary library file
unloadlibrary('libcantera_shared');
delete([cantera_root, ctName]);

disp('Cantera has been unloaded');
