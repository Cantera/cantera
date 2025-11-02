function mode = ctExecutionMode()
    % Retrieve the MATLAB toolbox execution mode
    if ~ctIsLoaded
        error("ctExecutionMode:NotLoaded", "Cantera toolbox is not loaded.")
    end

    hasGlobalCt = evalin('base','exist("ct","var") && isa(ct,"clibConfiguration")');
    if hasGlobalCt
        global ct
        mode = ct.ExecutionMode;
    else
        cfg = clibConfiguration("ctMatlab");
        mode = cfg.ExecutionMode;
    end
end
