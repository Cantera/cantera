function d = dataDirectories()
    % Get a cell array of the directories searched for data files. ::
    %
    %     >> ct.dataDirectories()
    %
    % Get a cell array of the directories Cantera searches for data files
    %
    % :return:
    %     Cell array with strings representing the data file search directories

    ct.isLoaded(true);
    d = ct.impl.getString('mCt_getDataDirectories', ';');
end
