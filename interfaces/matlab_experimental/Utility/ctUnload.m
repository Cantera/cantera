function ctUnload()
    % Unload the Cantear C Library from the Memory
    %
    if libisloaded(ctSharedLibrary)
        ctCleanUp;
        unloadlibrary(ctSharedLibrary);
    end

    disp('Cantera has been unloaded');
end
