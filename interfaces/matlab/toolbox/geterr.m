function e = geterr
try
    e = ctmethods(0,2); % getCanteraError;
catch ME
    e = getReport(ME);
end
