function output = getString(funcName, varargin)
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    persistent emptyCache stringCache maxLen
    if isempty(stringCache)
        stringCache = clib.array.ctMatlab.Char(256);
        maxLen = 256;
    end
    if isempty(emptyCache)
        emptyCache = clib.array.ctMatlab.Char(0);
    end

    buflen = clib.ctMatlab.(funcName)(varargin{:}, emptyCache);

    if buflen <= 0
        error('Cantera:ctError', ct.impl.getError());
    end

    % resize stringCache if needed
    if buflen > maxLen
        stringCache = clib.array.ctMatlab.Char(buflen);
        maxLen = buflen;
    end

    nchar = sum(cellfun(@ischar, varargin));
    if nchar == 0 || nchar == 1
        buflen = clib.ctMatlab.(funcName)(varargin{:}, stringCache);
    else
        error('not implemented - argument list contains more than one char array.')
    end

    buflen = double(buflen);
    ct.impl.checkErrorCode(buflen);

    % Discard the last character
    s = stringCache.double;
    s = s(1:buflen);
    if s(end) == 0
        s = s(1:end-1);
    end
    output = char(s);

end
