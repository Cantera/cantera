function output = ctString(funcName, varargin)
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    errorcode = [-1, -999.999, double(intmax('uint64'))];

    buf = clib.array.ctMatlab.Char(0);
    buflen = clib.ctMatlab.(funcName)(varargin{:}, buf);

    if buflen > 0
        buf = clib.array.ctMatlab.Char(buflen);

        nchar = sum(cellfun(@ischar, varargin));
        if nchar == 0 || nchar == 1
            iok = clib.ctMatlab.(funcName)(varargin{:}, buf);
        else
            error('not implemented - argument list contains more than one char array.')
        end
    else
        error('Cantera:ctError', ctGetErr);
    end

    iok = double(iok);
    if ismember(iok, errorcode)
        error('Cantera:ctError', ctGetErr);
    end

    output = char(buf.double);

end
