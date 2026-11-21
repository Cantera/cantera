function install(opts)
    % Set up the Cantera MATLAB toolbox after installing it. ::
    %
    %     >> ct.install
    %     >> ct.install(Force=true)
    %
    % Run once after installing the toolbox; afterwards use :func:`ct.load`.
    % This registers the data directory shipped with the toolbox and verifies
    % that Cantera loads. Settings persist as MATLAB preferences in the
    % "Cantera" group, and later calls do nothing unless the installation is
    % incomplete or ``Force`` is true. Use :func:`ct.uninstall` to clear them.
    %
    % :param DataDirectory:
    %     Data directory to register instead of the one shipped with the
    %     toolbox.
    % :param Force:
    %     Reinstall even if Cantera is already installed.

    arguments
        opts.DataDirectory (1,1) string = ""
        opts.Force (1,1) logical = false
    end

    GROUP = "Cantera";
    [interfaceDir, toolboxDir] = toolboxPaths();

    if ~opts.Force && isInstalled(GROUP, interfaceDir)
        fprintf("Cantera is already installed. " + ...
                "Use ct.install(Force=true) to reinstall.\n");
        return
    end

    dataDir = resolveDataDirectory(toolboxDir, opts.DataDirectory);
    setpref(GROUP, "DataDirectory", char(dataDir));
    fprintf("Registered Cantera data directory: %s\n", dataDir);

    % Set the flag only after verification, so a failed install is retried.
    verifyInstallation(dataDir);
    setpref(GROUP, "Installed", true);
    fprintf("Cantera %s installed successfully.\n", ct.version);
end

function tf = isInstalled(GROUP, interfaceDir)
    % The flag alone is not trusted: preferences outlive the toolbox when it is
    % removed in the Add-On Manager without running ct.uninstall.
    tf = ispref(GROUP, "Installed") && isequal(getpref(GROUP, "Installed"), true) ...
         && hasInterface(interfaceDir) && ispref(GROUP, "DataDirectory") ...
         && isfolder(string(getpref(GROUP, "DataDirectory")));
end

function dataDir = resolveDataDirectory(toolboxDir, override)
    if override ~= ""
        if ~isfolder(override)
            error("ct:install:BadDataDir", ...
                  "Data directory does not exist: %s", override);
        end
        dataDir = override;
        return
    end

    candidates = [
        fullfile(fileparts(toolboxDir), "data")             % packaged toolbox
        fullfile(fileparts(fileparts(toolboxDir)), "data")  % source checkout
    ];
    for c = candidates'
        if ~isempty(dir(fullfile(c, "*.yaml")))
            dataDir = c;
            return
        end
    end

    error("ct:install:NoDataDir", ...
          ("Could not find the Cantera data files shipped with the toolbox. " + ...
           "Pass a folder with ct.install(DataDirectory=<folder>)."));
end

function verifyInstallation(dataDir)
    if ct.isLoaded
        ct.addDataDirectories(dataDir);
    else
        ct.load();
    end

    registered = string(ct.dataDirectories());
    if ~any(arrayfun(@(d) samePath(d, dataDir), registered))
        error("ct:install:DataDirNotRegistered", ...
              "Cantera loaded, but did not register the data directory %s", ...
              dataDir);
    end
end

function tf = samePath(a, b)
    normalize = @(p) regexprep(strrep(string(p), "\", "/"), "/+$", "");
    if ispc
        tf = strcmpi(normalize(a), normalize(b));
    else
        tf = normalize(a) == normalize(b);
    end
end
