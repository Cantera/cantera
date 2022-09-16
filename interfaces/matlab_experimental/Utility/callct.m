function output = callct(varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    err1 = -1;
    err2 = -999.999;
    err3 = double(intmax('uint64'));

    funcName = varargin{1};
    output = calllib(ct, funcName, varargin{2:end});
    if output == err2 || output == err3
        output = err1;
    end
    if output == err1
        error(geterr);
    end
end
