function s = session()
    % Number of times the Cantera library has been unloaded in this MATLAB session.
    %
    % Handles issued by the library are only valid while this value is unchanged:
    % unloading frees every object the library holds, and a reloaded library reuses
    % handle values. The counter is stored on `groot` so that it survives
    % `clear all`, which would reset a persistent variable.

    s = getappdata(groot, 'CanteraSession');
    if isempty(s)
        s = 0;
    end
end
