function buildToolbox()
    % Get Cantera version from environment variable
    version = getenv('TOOLBOX_VERSION');
    if isempty(version)
        version = '0.0.0';  % fallback default
    end
    fprintf('Building toolbox version: %s\n', version);

    % Read the UUID for the toolbox
    uuidFile = fullfile(pwd, 'interfaces', 'matlab_experimental', 'Setup', ...
                        'Cantera_MATLAB_Toolbox.uuid');
    if isfile(uuidFile)
        guid = strtrim(fileread(uuidFile));
        fprintf('The unique identifier for the toolbox is: %s\n', guid);
    else
        error('A unique identifier for the toolbox does not exist!');
    end

    % Define output file
    outputFile = fullfile(pwd, 'interfaces', 'matlab_experimental', ...
                          ['Cantera_MATLAB_Toolbox_', version, '.mltbx']);

    % Create temporary folder for reorganizing and faster execution
    tmpDir = fullfile(pwd, 'tmpDir');
    mkdir(tmpDir);
    mapping = {
        'interfaces/matlab_experimental', 'toolbox';
        'samples/matlab_experimental', 'samples';
        'test/matlab_experimental', 'test/matlab_toolbox';
        'test/data', 'test/data';
        'data', 'data'
        };

    % Move all files to temporary folder
    allFiles = {};
    for i = 1:size(mapping, 1)
        src = fullfile(pwd, mapping{i, 1});
        dest = fullfile(tmpDir, mapping{i, 2});
        copyfile(src, dest);
        files = dir(fullfile(dest, '**', '*'));
        files = files(~[files.isdir]);
        filesList = fullfile({files.folder}, {files.name})';
        allFiles = [allFiles; filesList];
    end

    % Get relative paths
    allPaths = strsplit(genpath(tmpDir), pathsep);
    relPaths = cellfun(@(p) erase(p, [pwd filesep]), allPaths, 'UniformOutput', false);
    relPaths = relPaths(~cellfun(@isempty, relPaths));

    % Get path to the icon file
    iconFile = fullfile(pwd, 'doc', 'sphinx', '_static', 'images', 'cantera-logo.png');

    % Set up toolbox options
    opts = matlab.addons.toolbox.ToolboxOptions(tmpDir, guid);
    opts.ToolboxName                     = 'Cantera MATLAB Toolbox';
    opts.ToolboxVersion                  = version;
    opts.Summary                         = 'MATLAB interface for Cantera.';
    opts.Description                     = [
        'Cantera is an open-source suite of tools for problems involving', ...
        'involving chemical kinetics, thermodynamics, and transport processes.', ...
        'This toolbox includes the MATLAB interface for Cantera, example scripts, and data files.'
         ];
    opts.AuthorName                      = 'Cantera Developers'; % placeholder
    opts.AuthorEmail                     = 'developers@cantera.org'; % placeholder
    opts.ToolboxFiles                    = allFiles;
    opts.ToolboxMatlabPath               = relPaths;
    opts.ToolboxImageFile                = iconFile;
    opts.MinimumMatlabRelease            = 'R2014b';
    opts.OutputFile                      = outputFile;
    opts.SupportedPlatforms.Win64        = true;
    opts.SupportedPlatforms.Glnxa64      = true;
    opts.SupportedPlatforms.Maci64       = true;
    opts.SupportedPlatforms.MatlabOnline = false;
    % These options will be enabled when we host Cantera binaries somewhere
    % opts.RequiredAdditionalSoftware = [
    %     struct( ...
    %         "Name", "CanteraBinaries", ...
    %         "Platform", "win64", ...
    %         "DownloadURL", "placeholder for download url", ...
    %         "LicenseURL", "placeholder for license url"),
    %     struct( ...
    %         "Name", "CanteraBinaries", ...
    %         "Platform", "maci64", ...
    %         "DownloadURL", "placeholder for download url", ...
    %         "LicenseURL", "placeholder for license url"),
    %     struct( ...
    %         "Name", "CanteraBinaries", ...
    %         "Platform", "glnxa64", ...
    %         "DownloadURL", "placeholder for download url", ...
    %         "LicenseURL", "placeholder for license url"),
    %     ];

    % Package the toolbox
    try
        matlab.addons.toolbox.packageToolbox(opts);
        fprintf('✅ Toolbox built successfully!\n');
    catch ME
        fprintf('❌ Toolbox build failed: %s\n', ME.message);
    end

    % Remove the temporary folder
    rmdir(tmpDir, 's');
end
