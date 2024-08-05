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
    d = ctString('ct_getDataDirectories3', ';');
end
