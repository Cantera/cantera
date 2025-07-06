function downloadDependencies()
    % This file may not be necessary once we set up the MLTBX properly

    os = lower(computer);
    % Get platform-specific binaries
    switch os
        case {'pcwin64', 'pcwin'}
            platform = 'windows';
        case {'maci64'}
            platform = 'macos';
        case {'glnxa64', 'glnx86'}
            platform = 'linux';
        otherwise
            error('Unsupported platform: %s', os);
    end

    % Placeholder URL for ZIP files containing the binaries and headers
    baseUrl = 'https://cantera.com/downloads';
    zipName = sprintf('cantera-binaries-%s.zip', platform);
    url = fullfile(baseUrl, zipName);

    % Set target install folder in user path
    installDir = fullfile(userpath, 'Dependencies');
    if ~isfolder(installDir)
        fprintf('Downloading dependencies for %s...\n', platform);
        zipFile = fullfile(tempdir, zipName);
        websave(zipFile, url);
        unzip(zipFile, installDir);
        fprintf('Installed to: %s\n', installDir);
    else
        fprintf('Dependencies already installed at: %s\n', installDir);
    end

    % Add to path
    addpath(genpath(installDir));
end
