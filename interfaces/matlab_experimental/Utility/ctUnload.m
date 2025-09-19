function ctUnload()
    % Unload the Cantera C Library from the Memory

    global ct
    if ~isempty(ct)
        ctCleanUp;
        ct.unload
    end

    clear global ct

    disp('Cantera has been unloaded');
end
