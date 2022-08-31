function e = geterr()
    % GETERR
    % Get the error message from a Cantera error.
    %
    checklib;
    try
        buflen = calllib(ct, 'ct_getCanteraError', 0, '');
        aa = char(ones(1, buflen));
        ptr = libpointer('cstring', aa);
        [iok, bb] = calllib(ct, 'ct_getCanteraError', buflen, ptr);
        e = bb;
        clear aa bb ptr
    catch ME
        e = getReport(ME);
    end
end
