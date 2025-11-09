function e = getError()
    % Get the error message from a Cantera error.

    try
        e = ct.impl.getString('mCt_getCanteraError');
    catch ME
        e = getReport(ME);
    end

end
