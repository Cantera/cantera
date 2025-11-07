function ctUnload()
    % Unload Cantera MATLAB C++ interface and clean up MATLAB side
    if ~ctIsLoaded(false)
        return
    end

    try
        ctCleanUp;
    catch ME
        warning("ctUnload:CleanupFailed", ...
                "ctCleanUp failed (%s).", ME.message);
    end

    if ctExecutionMode() == "inprocess"
        warning("ctUnload:UnloadFailed", ...
                ("Unloading of `ctMatlab` library is not supported for " + ...
                 "'inprocess' execution mode. Restart MATLAB to unload."));
        return
    end

    global ct
    hasGlobalCt = evalin('base', ['exist("ct","var") && ', ...
                         'isa(ct,"matlab.cppclient.CLibraryConfiguration")']);
    if hasGlobalCt
        try
            ct.unload;
            disp("Cantera has been unloaded");
        catch ME
            warning("ctUnload:UnloadFailed", ...
                    "ct.unload failed (%s). Attempting fallback.", ME.message);
        end
        clear global ct
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
