function str = ctLib(libDir)
    % Return name of Cantera Shared Library depending on OS.
    function score = versionScore(libName)
        % Calculate sortable integer version score based on version number.
        % Pattern to match version numbers (x.y.z, x.y, or x)
        % Matches sequences of digits separated by dots
        pattern = '(\d+(?:\.\d+){0,2})(?=\.(?:dylib|so)|$)';

        % Find the version number
        tokens = regexp(libName, pattern, 'tokens', 'once');
        if ~isempty(tokens)
            verStr = tokens{1};
        else
            score = -1;
            return;
        end

        % Parse version string into comparable score, e.g. x00y00z
        parts = sscanf(verStr, '%d.%d.%d');
        parts(end+1:3) = 0;  % Pad to length 3
        score = parts(1)*1e6 + parts(2)*1e3 + parts(3);
    end

    if ispc
        if isempty(dir(fullfile(libDir, "cantera_shared.dll")))
            error('No matching shared libraries found in %s', libDir);
        end
        str = string(fullfile(libDir, "cantera_shared.dll"));
    elseif ismac
        % MATLAB requires exact library name `libcantera_shared.X.Y.Z.dylib` to
        % load correctly
        if ~isempty(dir(fullfile(libDir, 'libcantera_shared*.dylib')))
            baseName = 'libcantera_shared';
        elseif ~isempty(dir(fullfile(libDir, 'libcantera*.dylib')))
            baseName = 'libcantera';
        else
            error('No matching shared libraries found in %s', libDir);
        end
        files = dir(fullfile(libDir, [baseName, '*.dylib']));

        scores = arrayfun(@(f) versionScore(f.name), files);
        hasVersion = scores > 0;
        if any(hasVersion)
            [~, idx] = max(scores);
            str = string(fullfile(libDir, files(idx).name));
        else
            warning('No versioned library found, falling back to unversioned');
            str = string(fullfile(libDir, files(1).name));
        end
    elseif isunix
        % MATLAB requires exact library name `libcantera_shared.so.X` to
        % load correctly
        if ~isempty(dir(fullfile(libDir, 'libcantera_shared.so*')))
            baseName = 'libcantera_shared';
        elseif ~isempty(dir(fullfile(libDir, 'libcantera.so*')))
            baseName = 'libcantera';
        else
            error('No matching shared libraries found in %s', libDir);
        end
        files = dir(fullfile(libDir, [baseName, '*.so*']));

        scores = arrayfun(@(f) versionScore(f.name), files);
        hasVersion = scores > 0;
        if any(hasVersion)
            latestMajorVersion = scores >= round(max(scores) / 1e6) * 1e6;
            files = files(latestMajorVersion);
            scores = scores(latestMajorVersion);
            [~, idx] = min(scores);
            str = string(fullfile(libDir, files(idx).name));
        else
            warning('No versioned library found, falling back to unversioned');
            str = string(fullfile(libDir, files(1).name));
        end
    end
end
