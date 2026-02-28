function unload()
    % Unload Cantera MATLAB C++ interface and clean up MATLAB side
    if ~ct.isLoaded(false)
        return
    end

    try
        ct.cleanUp;
    catch ME
        warning("ct.unload:CleanupFailed", ...
                "cleanUp failed (%s).", ME.message);
    end

    % Clear the persistent cache even if cleanUp failed
    clear getString getArray

    if ct.executionMode() == "inprocess"
        warning("ct.unload:UnloadFailed", ...
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
            warning("ct.unload:UnloadFailed", ...
                    "ctMatlab.unload failed (%s). Attempting fallback.", ME.message);
        end
        clear global ctMatlab
    else
        cfg = clibConfiguration("ctMatlab");
        try
            unload(cfg);
            disp("Cantera has been unloaded");
        catch ME
            warning("ct.unload:UnloadFailed", ...
                "unload(clibConfiguration) failed (%s).", ME.message);
        end
    end
end
