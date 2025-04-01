% Read the file contents
fin = 'definectMatlabInterface.m'
fout = 'tmp.m'
fid = fopen(fin, 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = lines{1};

% Process each line to find and replace <SHAPE>
for i = 1:length(lines)
    if contains(lines{i}, '%defineArgument') && contains(lines{i}, '<SHAPE>')
        % Extract argument names from previous defineArgument lines
        j = i - 1;
        while j > 0 && ~contains(lines{j}, 'defineArgument')
            j = j - 1;
        end
        
        if j > 0
            % Extract previous argument name
            tokens = regexp(lines{j}, 'defineArgument\(.*?,\s*"(.*?)"', 'tokens');
            if ~isempty(tokens)
                shapeArg = tokens{1}{1};
                % Replace <SHAPE> with "shapeArg" in the current line, preserving comments
                lines{i} = regexprep(lines{i}, '<SHAPE>', ['"' shapeArg '"']);
            end
        end
        
        % Find the beginning of the function block (%% marker)
        k = i;
        while k > 0 && ~startsWith(lines{k}, '%%')
            k = k - 1;
        end

        l = i + 1;
        while l <= length(lines) && ~startsWith(lines{l}, '%%')
            l = l + 1;
        end
        
        % Uncomment all lines except the first two in the block, preserving comments
        for m = k+2:l-1
            if startsWith(strtrim(lines{m}), '%')
                comment_idx = strfind(lines{m}, '%');
                if ~isempty(comment_idx)
                    lines{m} = [erase(lines{m}(1:comment_idx(1)), '%'), lines{m}(comment_idx(1) + 1:end)];
                end
            end
        end
    end
end

% Write the modified content back
fid = fopen(fout, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
