function phases = ImportPhases(src, phasenames)
    % ImportPhases ::
    %
    % phases = ImportPhases(src, phasenames)
    %
    % :param src:
    %      YAML file containing the interface or edge phase.
    % :param phasenames:
    %      Name of a single phase (char/string) or a cell array of phase names.
    % :return:
    %      A cell array of :mat:class:`Solution` objects created from each phase.

    isLoaded(true);

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
