function d = dataDirectories()
    % Get a cell array of the directories searched for data files. ::
    %
    %     >> ct.dataDirectories()
    %
    % Get a cell array of the directories Cantera searches for data files
    %
    % :return:
    %     Cell array with strings representing the data file search directories

    isLoaded(true);
    d = ctString('mCt_getDataDirectories', ';');
end
