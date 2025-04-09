function i = ctIsLoaded()
    i = false;

    global ct
    if ~isempty(ct)
        i = ct.Loaded;
    end

    if ~i
        error('Cantera is not loaded.');
    end
end
