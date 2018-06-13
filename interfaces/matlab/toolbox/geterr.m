function e = geterr
% GETERR  Get the error message from a Cantera error.
% e = geterr
%

try
    e = ctmethods(0,2); % getCanteraError;
catch ME
    e = getReport(ME);
end
