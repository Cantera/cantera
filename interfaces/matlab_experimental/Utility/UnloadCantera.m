% Unload the Cantear C Library from the Memory
%
if libisloaded(ct)
    unloadlibrary(ct);
end

disp('Cantera has been unloaded');
