function str = ct()
    % Return name of 'cantera_shared' library.
    if ispc
        str = 'cantera_shared';
    else
        str = 'libcantera_shared';
    end

end
