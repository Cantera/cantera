function ctEditLibraryDefinitions(fileDir)
    fname = fileDir + "/definectMatlab.m";
    % Read original file
    lines = string(readlines(fname));

    % Find block boundaries
    blockStartIdx = find(startsWith(lines, "%% "));
    blockEndIdx = [blockStartIdx(2:end) - 1; numel(lines)];

    % Process each block
    for i = 2:numel(blockStartIdx)
        blockLines = lines(blockStartIdx(i):blockEndIdx(i));

        if any(contains(blockLines, '<DIRECTION>'))
            continue
        end

        if any(contains(blockLines, '<SHAPE>'))
            for j = 3:numel(blockLines)
                if contains(blockLines(j), '<SHAPE>')
                    k = j - 1;
                    shape = '';

                    if (~contains(blockLines(k-1), 'defineArgument') && ...
                        ~contains(blockLines(j), '.Char') && ...
                        contains(blockLines(k), 'int32')) || ...
                        contains(blockLines(k), 'clib.array') || ...
                        contains(blockLines(k), 'double')
                        % Resolve edge cases
                        shape = '2';
                    elseif any(contains(blockLines, 'trans_getMultiDiffCoeffsDefinition')) || ...
                           any(contains(blockLines, 'trans_getBinDiffCoeffsDefinition')) && ...
                           contains(blockLines(j), '"d"')
                            shape = '["ld","ld"]';
                    else
                        % Use name of preceding scalar variable
                        tokens = regexp(blockLines(k), '"(\w+)"', 'tokens', 'once');
                        if ~isempty(tokens)
                            shape = ['"' tokens{1} '"'];
                        else
                            shape = '1';
                        end
                    end

                    % Replace <SHAPE> with determined value
                    blockLines(j) = replace(blockLines(j), '<SHAPE>', shape);
                end
            end
        end

        % Uncomment and clean up trailing comments
        for j = 3:numel(blockLines)
            line = strtrim(blockLines(j));
            if startsWith(line, "%")
                line = extractAfter(line, 1);
            end
            pctIdx = strfind(line, '%');
            if ~isempty(pctIdx)
                lastPct = pctIdx(end);
                line = strtrim(extractBefore(line, lastPct));
            end
            blockLines(j) = line;
        end

        lines(blockStartIdx(i):blockEndIdx(i)) = blockLines;
    end

    % Write to libDef file
    writelines(lines, fname);
end
