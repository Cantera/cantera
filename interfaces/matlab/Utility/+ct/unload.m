function unload()
    % Unload Cantera MATLAB C++ interface and clean up MATLAB side
    if ~isLoaded(false)
        return
    end

    try
        cleanUp;
    catch ME
        warning("ctUnload:CleanupFailed", ...
                "cleanUp failed (%s).", ME.message);
    end

    if executionMode() == "inprocess"
        warning("ctUnload:UnloadFailed", ...
                ("Unloading of `ctMatlab` library is not supported for " + ...
                 "'inprocess' execution mode. Restart MATLAB to unload."));
        return
    end

    global ctMatlab
    hasGlobalCt = evalin('base', ['exist("ctMatlab","var") && ', ...
                         'isa(ctMatlab,"matlab.cppclient.CLibraryConfiguration")']);
    if hasGlobalCt
        try
            ctMatlab.unload;
            disp("Cantera has been unloaded");
        catch ME
            warning("ctUnload:UnloadFailed", ...
                    "ctMatlab.unload failed (%s). Attempting fallback.", ME.message);
        end
        clear global ctMatlab
    else
        cfg = clibConfiguration("ctMatlab");
        try
            unload(cfg);
            disp("Cantera has been unloaded");
        catch ME
            warning("ctUnload:UnloadFailed", ...
                "unload(clibConfiguration) failed (%s).", ME.message);
        end
    end
end
