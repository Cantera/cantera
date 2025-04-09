% Load the MATLAB-to-C++ Interface in Out-of-process Mode
global ct
if isempty(ct)
    ct = clibConfiguration("ctMatlab", ExecutionMode = "outofprocess");
end

disp(sprintf('Cantera %s is ready for use.', ctVersion));
