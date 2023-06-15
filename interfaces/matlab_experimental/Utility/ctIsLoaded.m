function i = ctIsLoaded()
    if ~libisloaded(ctLib)
        error('Cantera is not loaded');
    end

    i = true;

end
