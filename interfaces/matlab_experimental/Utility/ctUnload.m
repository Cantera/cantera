function ctUnload()
    %UNLOAD Cantera C++ interface and clean up MATLAB side

    try
        ctCleanUp;
    catch ME
        warning("ctUnload:CleanupFailed", ...
                "ctCleanUp failed (%s).", ME.message);
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
