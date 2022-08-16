function d = getDataDirectories()
    % Get a cell array of the directories searched for data files.
    % getdatadirectories()
    % Get a cell array of the directories Cantera searches for data files
    %
    % :return:
    %     Cell array with strings representing the data file search directories
    %
    checklib;
    buflen = callct('ct_getDataDirectories', 0, '', ';');
    aa = char(zeros(1, buflen));
    [~, aa, ~] = callct('ct_getDataDirectories', buflen, aa, ';');
    d = aa;
end
