function i = ctIsLoaded()
    i = false;

    hasGlobalCt = evalin('base','exist("ct","var") && isa(ct,"clibConfiguration")');
    if hasGlobalCt
        global ct
        i = ct.Loaded;
    else
        cfg = clibConfiguration("ctMatlab");
        i = cfg.Loaded;
    end

    if ~i
        error('Cantera is not loaded.');
    end
end
