function output = callct(varargin)
    % CALLCT
    % This is a simplified single output variant
    err1 = -1;
    err2 = -999.999;
    err3 = double(intmax('uint64'));

    methodname = varargin{1};
    output = calllib(ct, methodname, varargin{2:end});
    if output == err2 || output == err3
        output = -1;
    end
    if output == -1
        geterr;
    end
end
