function tf = hasInterface(folder)
    % True if folder contains the compiled interface for this platform.
    %
    % The packaged toolbox ships every platform's interface in one folder, so
    % only the current platform's file counts.
    arguments
        folder (1,1) string
    end

    if ispc
        ext = "dll";
    elseif ismac
        ext = "dylib";
    else
        ext = "so";
    end
    tf = isfile(fullfile(folder, "ctMatlabInterface." + ext));
end
