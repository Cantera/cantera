function uninstall()
    % Remove all persisted Cantera settings. ::
    %
    %     >> ct.uninstall
    %
    % Clears the "Cantera" MATLAB preference group written by :func:`ct.install`.
    % The toolbox itself is removed in the Add-On Manager; run this first, since
    % this function is part of the toolbox.

    % A loaded library keeps its files locked on Windows, which would stop the
    % Add-On Manager from deleting them.
    if ct.isLoaded
        ct.unload();
    end

    if ispref("Cantera")
        rmpref("Cantera");
        fprintf("Removed Cantera MATLAB preferences.\n");
    else
        fprintf("No Cantera MATLAB preferences to remove.\n");
    end
end
