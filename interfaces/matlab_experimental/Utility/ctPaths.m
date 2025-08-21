function paths = ctPaths(varargin)
    % ctPaths ::
    %   Manage Cantera MATLAB toolbox paths
    %
    %   >> ctPaths()
    %       Get current config as struct
    %   >> ctPaths('Set', envRoot)
    %       Configure paths based on the conda environment where Cantera is installed
    %   >> ctPaths('Unset')
    %       Unset saved paths
    %
    % :return:
    %   paths: struct with fields 'libPath', 'includePath', and
    %   'toolboxPath'

    prefGroup = 'Cantera';
    prefName = 'Paths';
    if nargin == 0
        if ispref(prefGroup, prefName)
            paths = getpref(prefGroup, prefName);
        else
            paths = struct();
        end
        return
    end

    option = validatestring(varargin{1}, {'Set', 'Unset'}, mfilename, 'option', 1);

    switch option
        case 'Set'
            if numel(varargin) < 2
                error('ctPaths:SetMissing', ...
                      'You must provide the conda environment root with ctPaths(''Set'', envRoot).');
            end
            envRoot = varargin{2};

            if ispc
                libPath     = fullfile(envRoot, 'Library', 'bin');
                includePath = fullfile(envRoot, 'Library', 'include');
            else
                libPath     = fullfile(envRoot, 'lib');
                includePath = fullfile(envRoot, 'include');
            end

            samplePath = fullfile(envRoot, 'share', 'cantera', 'samples');
            dataPath   = fullfile(envRoot, 'share', 'cantera', 'data');

            thisFile    = mfilename('fullpath');
            toolboxPath = fileparts(fileparts(thisFile));

            paths = struct( ...
                'libPath',     libPath, ...
                'includePath', includePath, ...
                'samplePath',  samplePath, ...
                'dataPath',    dataPath, ...
                'toolboxPath', toolboxPath );

            setpref(prefGroup, prefName, paths);

            addpath(genpath(toolboxPath), genpath(samplePath));
            savepath;

        case 'Unset'
            if ispref(prefGroup, prefName)
                savedPaths = getpref(prefGroup, prefName);
                rmpref(prefGroup, prefName);
            else
                paths = struct();
            end

            if isfield(savedPaths, 'toolboxPath') && isfolder(savedPaths.toolboxPath)
                rmpath(genpath(savedPaths.toolboxPath));
            end

            if isfield(savedPaths, 'samplePath') && isfolder(savedPaths.samplePath)
                rmpath(genpath(savedPaths.samplePath));
            end

            savepath;
    end
end
