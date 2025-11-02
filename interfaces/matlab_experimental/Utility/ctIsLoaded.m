function i = ctIsLoaded(throw)
    % Check whether Cantera MATLAB C++ interface is loaded.
    %
    %     >> s = ctIsLoaded
    %     >> s = ctIsLoaded(true)
    %
    % :param throw:
    %     If interface is not loaded, throw error instead of issuing warning;
    %     default is ``false``.
    % :return:
    %     ``true`` if loaded, otherwise ``false``.

    arguments
        throw (1,1) logical = false
    end
    i = false;

    hasGlobalCt = evalin('base', 'exist("ct","var") && isa(ct,"clibConfiguration")');
    if hasGlobalCt
        global ct
        i = ct.Loaded;
    else
        cfg = clibConfiguration("ctMatlab");
        i = cfg.Loaded;
    end

    if ~i && throw
        error('Cantera is not loaded.');
    end
end
