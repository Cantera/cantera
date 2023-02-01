function str = ctLib()
    % Return name of Cantera Shared Library depending on OS.
    if ispc
        str = 'cantera_shared';
    else
        str = 'libcantera_shared';
    end

end
