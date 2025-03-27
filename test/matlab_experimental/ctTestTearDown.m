clear all
% Unload Cantera
if ispc
    lib = 'cantera_shared';
else
    lib = 'libcantera_shared';
end

unloadlibrary(lib);
disp('Cantera has been unloaded');
