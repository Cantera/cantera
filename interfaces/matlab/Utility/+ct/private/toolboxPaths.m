function [interfaceDir, toolboxDir] = toolboxPaths()
    % Locate the compiled interface folder and the toolbox root.
    %
    % Paths are derived from this file's location, so they hold both for the
    % packaged toolbox and for a source checkout (interfaces/matlab).
    %
    % :return interfaceDir:
    %     Folder holding the compiled interface and the Cantera shared library.
    % :return toolboxDir:
    %     Toolbox root containing the +ct package and the Utility folder.

    ctDir = fileparts(fileparts(mfilename("fullpath")));
    toolboxDir = string(fileparts(fileparts(ctDir)));
    interfaceDir = fullfile(toolboxDir, "+ct", "+impl", "ctMatlab");
end
