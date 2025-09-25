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

            thisFile    = mfilename('fullpath');
            toolboxPath = fileparts(fileparts(thisFile));

            paths = struct( ...
                'libPath',     libPath, ...
                'includePath', includePath, ...
                'toolboxPath', toolboxPath );

            setpref(prefGroup, prefName, paths);

            addpath(genpath(toolboxPath));

        case 'Unset'
            if ispref(prefGroup, prefName)
                rmpref(prefGroup, prefName);
            else
                paths = struct();
            end

    end
end
