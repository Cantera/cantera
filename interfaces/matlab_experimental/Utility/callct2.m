function output = callct2(varargin)
    % CALLCT2
    % Calls Cantera library functions with string outputs and returns
    % errors if necessary.

    err1 = -1;

    buflen = calllib(ct, varargin, 0, '');
    if buflen > 0
        aa = char(ones(1, buflen));
        ptr = libpointer('cstring', aa);
        [iok, bb] = calllib(ct, varargin, buflen, ptr);
        output = bb;
        clear aa bb ptr;
    else
        error(geterr);
    end
    if iok == -err1
        error(geterr);
    end
end
