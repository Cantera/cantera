% get list of directories
% dirs = strread(path, '%s', 'delimiter', pathsep);
dirs = regexp(path, ['([^' pathsep ']*)'], 'match');

% if 'cantera' is already in the path, we want to remove it
for i = 1:length(dirs)
    if strfind(dirs{i}, 'Cantera')
        rmpath(dirs{i});
        continue;
    end
    if strfind(dirs{i}, 'cantera')
        rmpath(dirs{i});
    end
end

path(path, [pwd filesep 'toolbox']);
path(path, [pwd filesep 'toolbox' filesep '1D']);

% Set path to Python module
%setenv('PYTHONPATH', [pwd filesep 'interfaces' filesep 'python'])
