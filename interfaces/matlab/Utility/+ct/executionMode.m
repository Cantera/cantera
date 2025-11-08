function mode = executionMode()
    % Retrieve the MATLAB toolbox execution mode
    if ~ct.isLoaded
        error("executionMode:NotLoaded", "Cantera toolbox is not loaded.")
    end

    hasGlobalCt = evalin('base','exist("ctMatlab","var") && isa(ctMatlab,"clibConfiguration")');
    if hasGlobalCt
        global ctMatlab
        mode = ctMatlab.ExecutionMode;
    else
        cfg = clibConfiguration("ctMatlab");
        mode = cfg.ExecutionMode;
    end
end
