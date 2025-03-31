function output = ctString(funcName, varargin)
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    err1 = -1;

    [buflen, ~] = clib.ctMatlab.(funcName)(varargin{:}, 0);

    if buflen > 0
        % Convert all buflen to double for better data type compatibility
        buflen = double(buflen);

        nchar = sum(cellfun(@ischar, varargin));
        if nchar == 0 || nchar == 1
            [iok, str] = clib.ctMatlab.(funcName)(varargin{:}, buflen);
        else
            error('not implemented - argument list contains more than one char array.')
        end
    else
        error('Cantera:ctError', ctGetErr);
    end

    iok = double(iok);
    if iok == err1
        error('Cantera:ctError', ctGetErr);
    end

    output = str;

end
