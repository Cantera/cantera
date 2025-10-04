function str = ctLib(libDir)
    % Return name of Cantera Shared Library depending on OS.
    if ispc
        str = libDir + "/cantera_shared.dll";
    elseif ismac
        if ~isempty(dir(fullfile(libDir, 'libcantera.dylib')))
            baseName = 'libcantera';
        elseif ~isempty(dir(fullfile(libDir, 'libcantera_shared.dylib')))
            baseName = 'libcantera_shared';
        else
            error('No matching shared libraries found in %s', libDir);
        end
        files = dir(fullfile(libDir, [baseName, '*.dylib']));
        versioned = {};
        for k = 1:numel(files)
            name = files(k).name;
            tokens = regexp(name, [baseName,'\.([0-9]+(\.[0-9]+)*)\.dylib$'], 'tokens', 'once');
            if ~isempty(tokens)
                versioned{end+1} = name;
            end
        end
        if isempty(versioned)
            warning('No versioned library found, falling back to unversioned');
            str = libDir + baseName + ".dylib";
        else
            [~, idx] = max(cellfun(@(s) sscanf(s, [baseName, '.%f.dylib']), versioned));
            str = libDir + "/" + versioned{idx};
        end
    elseif isunix
        if ~isempty(dir(fullfile(libDir, 'libcantera.so*')))
            baseName = 'libcantera';
        elseif ~isempty(dir(fullfile(libDir, 'libcantera_shared.so*')))
            baseName = 'libcantera_shared';
        else
            error('No matching shared libraries found in %s', libDir);
        end
        files = dir(fullfile(libDir, [baseName, '.so*']));
        versioned = {};
        for k = 1:numel(files)
            name = files(k).name;
            tokens = regexp(name, '\.so\.(\d+)(\.\d+)*$', 'once');
            if ~isempty(tokens)
                versioned{end+1} = name;
            end
        end
        if isempty(versioned)
            warning('No versioned library found, falling back to unversioned');
            str = libDir + baseName + ".so";
        else
            [~, idx] = min(cellfun(@(s) sscanf(s, [baseName, '.so.%f']), versioned));
            str = libDir + "/" + versioned{idx};
        end
    end

end
