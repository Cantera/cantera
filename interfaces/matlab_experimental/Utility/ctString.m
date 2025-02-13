function output = ctString(funcName, varargin)
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    err1 = -1;

    buflen = calllib(ctLib, funcName, varargin{:}, 0, '') + 1;

    if buflen > 0
        aa = char(ones(1, buflen));
        ptr = libpointer('cstring', aa);
        nchar = sum(cellfun(@ischar, varargin));
        if nchar == 0
            % varargin does not contain char array
            [iok, bb] = calllib(ctLib, funcName, varargin{:}, buflen, ptr);
        elseif nchar == 1
            % varargin contains one char array, which MATLAB returns in second place
            [iok, ~, bb] = calllib(ctLib, funcName, varargin{:}, buflen, ptr);
        else
            error('not implemented - argument list contains more than one char array.')
        end
        output = bb;
        clear aa bb ptr;
    else
        error('Cantera:ctError', ctGetErr);
    end

    if iok == err1
        error('Cantera:ctError', ctGetErr);
    end

end
