function ctUnload()
    % Unload the Cantera C Library from the Memory

    if libisloaded(ctLib)
        ctCleanUp;
        unloadlibrary(ctLib);
    end

    disp('Cantera has been unloaded');
end
