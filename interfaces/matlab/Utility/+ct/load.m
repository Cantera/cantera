function load(mode)
    % Load the MATLAB-to-C++ Interface
    %
    % Usage:
    %   ct.load                % defaults to 'outofprocess'
    %   ct.load('inprocess')   % load in-process
    %   ct.load('outofprocess')% load out-of-process
    %
    % After installing the toolbox, run :func:`ct.install` once before loading.

    arguments
        mode (1,1) string {mustBeMember(mode, ...
                           ["inprocess", "outofprocess"])} = "outofprocess"
    end

    if ct.isLoaded
        mode_ = ct.executionMode;
        msg = sprintf('Cantera %s is already loaded (%s mode).', ct.version, mode_);
        warning("load:IsLoaded", msg);
        if mode == mode_ && mode == "inprocess"
            return
        elseif mode ~= mode_ && mode_ == "inprocess"
            error("load:LoadFailed", ...
                  ("Unloading of `ctMatlab` library is not supported for " + ...
                   "'inprocess' execution mode. Restart MATLAB to update."));
            return
        else
            ct.unload();
        end
    end

    libDir = toolboxPaths();
    if ~hasInterface(libDir)
        error("ct:load:NotInstalled", ...
              ("The Cantera MATLAB interface for this platform was not found " + ...
               "in %s.\nBuild it with ct.buildInterface, or reinstall the " + ...
               "toolbox."), libDir);
    end
    if ~any(strcmp(strsplit(path, pathsep), libDir))
        addpath(libDir);
    end

    % The interface folder also holds the Cantera shared library, which the OS
    % loader must find. Out-of-process mode on R2025a+ can pass environment
    % variables to the child library process; otherwise set them here.
    useOOPEnv = ~isMATLABReleaseOlderThan("R2025a") && mode == "outofprocess";
    pathVar = buildLoaderPathVar(libDir);
    if ~useOOPEnv
        applyLoaderEnvInProcess(pathVar);
    end

    global ctMatlab
    if ~ct.isLoaded
        if useOOPEnv && numEntries(pathVar) > 0
            ctMatlab = clibConfiguration("ctMatlab", ExecutionMode=mode, ...
                                   OutOfProcessEnvironmentVariables=pathVar);
        else
            ctMatlab = clibConfiguration("ctMatlab", ExecutionMode=mode);
        end
    end

    ctVersion = ct.version;
    addInstalledDataDirectory();

    fprintf('Cantera %s is ready for use (%s mode).\n', ctVersion, mode);
end

function pathVar = buildLoaderPathVar(libDir)
    % Build the OS-specific dynamic-loader environment dictionary, prepending
    % the library directory (and, on macOS, MATLAB's own libraries) to the
    % current value of the relevant environment variable.
    %
    % Nothing is set on Linux: the interface finds the Cantera library in its
    % own folder through its RUNPATH ($ORIGIN). Passing MATLAB's
    % LD_LIBRARY_PATH on to the library process would make Cantera's
    % dependencies resolve to the copies bundled with MATLAB instead.
    pathVar = dictionary(string.empty, string.empty);

    if ispc
        varName = "PATH";
        extra = strings(0, 1);
    elseif ismac
        varName = "DYLD_LIBRARY_PATH";
        arch = computer("arch");
        extra = [matlabroot + "/bin/" + arch; matlabroot + "/sys/os/" + arch];
    else
        return
    end

    parts = [libDir; extra];
    existing = string(getenv(varName));
    if existing ~= ""
        parts = [parts; split(existing, pathsep)];
    end
    parts = unique(parts(parts ~= ""), "stable");

    pathVar(varName) = strjoin(parts, pathsep);
end

function applyLoaderEnvInProcess(pathVar)
    % Apply the loader environment via setenv. Used when
    % OutOfProcessEnvironmentVariables is unavailable (in-process mode, or
    % MATLAB older than R2025a).
    k = keys(pathVar);
    if isempty(k)
        return
    end
    for i = 1:numel(k)
        setenv(char(k(i)), char(pathVar(k(i))));
    end

    if ~ispc
        warning("ct:load:InProcessEnvUnreliable", ...
                ("Setting %s after MATLAB has started may not affect the " + ...
                 "dynamic loader. If Cantera fails to load, use out-of-process " + ...
                 "mode on R2025a or newer, or launch MATLAB with the library " + ...
                 "directory already on the loader path."), char(k(1)));
    end
end

function addInstalledDataDirectory()
    % Data directories live in the loaded library, so the one registered by
    % ct.install has to be re-added on every load.
    GROUP = "Cantera";
    if ~ispref(GROUP, "DataDirectory")
        return
    end

    dataDir = string(getpref(GROUP, "DataDirectory"));
    dataDir = reshape(dataDir, 1, []);
    missing = dataDir(~arrayfun(@isfolder, dataDir));
    if ~isempty(missing)
        warning("ct:load:MissingDataDir", ...
                ("Cantera data directory no longer exists: %s\n" + ...
                 "Run ct.install(Force=true) to reconfigure."), ...
                strjoin(missing, ", "));
    end

    dataDir = dataDir(arrayfun(@isfolder, dataDir));
    if ~isempty(dataDir)
        ct.addDataDirectories(dataDir);
    end
end
