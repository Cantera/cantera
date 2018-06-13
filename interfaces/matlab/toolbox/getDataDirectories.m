function d = getDataDirectories()
% GETDATADIRECTORIES  Get a cell array of the directories searched for data files.
% getdatadirectories()
% Get a cell array of the directories Cantera searches for data files
%
% :return:
%     Cell array with strings representing the data file search directories
%

d = strsplit(ctmethods(0, 5, ';'), ';');
