function paths = ctPaths(varargin)
    % ctPaths ::
    % Configure or retrieve the library/header/toolbox paths for Cantera.
    % The paths are stored as MATLAB preferences.
    %
    %   >> ctPaths()                       % Get current config as struct
    %   >> ctPaths(libPath, includePath)   % Set library/header/toolbox paths
    %   >> ctPaths('clear')                % Clear saved paths
    %
    % :return:
    %   paths: struct with fields 'libPath', 'includePath', and
    %   'toolboxPath'

    if nargin == 1 && strcmp(varargin{1}, 'clear')
        if ispref('Cantera', 'Paths')
            paths = getpref('Cantera', 'Paths');
            subDirs = strsplit(genpath(paths.toolboxPath), pathsep);
            currentPaths = strsplit(path, pathsep);
            subdirsToRemove = intersect(subDirs, currentPaths);
            if ~isempty(subdirsToRemove)
                rmpath(subdirsToRemove{:});
            end
            rmpref('Cantera', 'Paths');
        end
        paths = struct('libPath', '', 'includePath', '', 'toolboxPath', '');
        return
    elseif nargin == 3 && all(cellfun(@ischar, varargin))
        paths = struct('libPath', varargin{1}, ...
                       'includePath', varargin{2}, ...
                       'toolboxPath', varargin{3});
        setpref('Cantera', 'Paths', paths);
    else
        % Load from saved preferences if available
        if ispref('Cantera', 'Paths')
            paths = getpref('Cantera', 'Paths');
            return
        else
            paths = struct('libPath', '', 'includePath', '', 'toolboxPath', '');
        end
    end

    if any(cellfun(@isempty, {paths.libPath, paths.includePath, paths.toolboxPath}))
        error('ctPaths:MissingPath', ...
              'Library, header, and toolbox paths must be configured with ctPaths(libPath,includePath, toolboxPath).');
    end

    mapping = {
        'interfaces/matlab_experimental', 'toolbox';
        'samples/matlab_experimental', 'samples';
        'test/matlab_experimental', 'test/matlab_toolbox';
        'test/data', 'test/data';
        'data', 'data'
        };

    % Check whether user is using Cantera source code or MLTBX based on folder structure
    if isfolder(fullfile(paths.toolboxPath, 'interfaces'))
        col = 1;
    else
        col = 2;
    end

    for i = 1:size(mapping, 1)
        subdir = fullfile(paths.toolboxPath, mapping{i, col});
        if isfolder(subdir)
            addpath(genpath(subdir));
        else
            warning('ctPaths:MissingDirectory', ...
                    'Directory not found: %s', subdir);
        end
    end

    savepath();
end
