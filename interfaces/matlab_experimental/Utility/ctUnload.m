function ctUnload()
    % Unload the Cantear C Library from the Memory
    %
    if libisloaded(ctLib)
        ctCleanUp;
        unloadlibrary(ctLib);
    end

    disp('Cantera has been unloaded');
end
