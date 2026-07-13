function addDataDirectories(dirs)
    % Add one or more directories to the data file search path. ::
    %
    %     >> ct.addDataDirectories('/path/to/data')
    %     >> ct.addDataDirectories(["/path/one", "/path/two"])
    %
    % Directories are added to the front of the Cantera data file search path,
    % so the most recently added directory is searched first. Use
    % :mat:func:`dataDirectories` to inspect the current search path.
    %
    % :param dirs:
    %     A string, char array, or string array of directories to add to the
    %     Cantera data file search path.
    arguments
        dirs (1,:) string
    end

    ct.isLoaded(true);
    for d = dirs
        ct.impl.call('mCt_addDataDirectory', char(d));
    end
end
