function e = ctGetErr()
    % Get the error message from a Cantera error.

    try
        e = ctString('mCt_getCanteraError');
    catch ME
        e = getReport(ME);
    end

end
