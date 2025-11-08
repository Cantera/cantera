function i = isLoaded(throw)
    % Check whether Cantera MATLAB C++ interface is loaded.
    %
    %     >> s = ct.isLoaded
    %     >> s = ct.isLoaded(true)
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

    hasGlobalCt = evalin('base', 'exist("ctMatlab","var") && isa(ctMatlab,"clibConfiguration")');
    if hasGlobalCt
        global ctMatlab
        i = ctMatlab.Loaded;
    else
        cfg = clibConfiguration("ctMatlab");
        i = cfg.Loaded;
    end

    if ~i && throw
        error('Cantera is not loaded.');
    end
end
