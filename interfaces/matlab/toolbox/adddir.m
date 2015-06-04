function adddir(d)
% ADDDIR  Add a directory to the search path.
% adddir(d)
% Adds directory ``d`` to the set of directories where Cantera looks for
% input and data files.
%
% :param d:
%     Path to add to the MATLAB search path.
%

ctmethods(0, 3, d);
