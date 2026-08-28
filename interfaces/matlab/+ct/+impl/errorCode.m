function err = errorCode()
    % Values returned by CLib functions to signal that an exception was raised.
    %
    % `ERR` (-999) and `DERR` (-999.999) are the sentinels defined by
    % `clib_defs.h`; -1 and -2 are returned in their place by generated functions
    % whose normal return value is a length or a handle. `intmax('uint64')` is the
    % C++ `npos`.

    err = [-1, -2, -999, -999.999, double(intmax('uint64'))];
end
