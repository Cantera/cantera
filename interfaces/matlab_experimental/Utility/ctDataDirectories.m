function d = ctDataDirectories()
    % Get a cell array of the directories searched for data files. ::
    %
    %     >> ctDataDirectories()
    %
    % Get a cell array of the directories Cantera searches for data files
    %
    % :return:
    %     Cell array with strings representing the data file search directories

    ctIsLoaded;
    buflen = calllib(ct, 'ct_getDataDirectories', 0, '', ';');
    aa = char(ones(1, buflen));
    ptr = libpointer('cstring', aa);
    [~, bb, ~] = calllib(ct, 'ct_getDataDirectories', buflen, ptr, ';');
    d = bb;
    clear aa, bb, ptr
end
