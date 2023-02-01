function output = ctString(varargin)
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    err1 = -1;

    funcName = varargin{1};
    buflen = calllib(ctLib, funcName, varargin{2:end}, 0, '');

    if buflen > 0
        aa = char(ones(1, buflen));
        ptr = libpointer('cstring', aa);
        [iok, bb] = calllib(ctLib, funcName, varargin{2:end}, buflen, ptr);
        output = bb;
        clear aa bb ptr;
    else
        error(ctGetErr);
    end

    if iok == -err1
        error(ctGetErr);
    end

end
