function phases = ctImportPhases(src, phasenames)
    ctIsLoaded;

    if nargin < 2
        error('Please specify the source file and list of phases to import');
    end

    phases = {};

    if ischar(phasenames)
        phases = {Solution(src, phasenames)};
        return
    elseif iscell(phasenames)
        for name = phasenames
            phases = [phases; {Solution(src, name{:})}];
        end
    else
        error (['Invalid type for phasenames, ' , ...
               'expected string or a cell array of strings']);
    end

end
