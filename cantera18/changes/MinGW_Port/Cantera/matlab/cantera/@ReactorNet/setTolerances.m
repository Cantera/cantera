function setTolerances(r, rtol, atol)
% SETTOLERANCES - Set error tolerances.
%
% rtol - scalar relative error tolerance
% atol - scalar absolute error tolerance
%
reactornetmethods(7, reactornet_hndl(r), rtol, atol);
