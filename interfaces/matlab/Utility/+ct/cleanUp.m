function cleanUp()
    % Delete all stored Cantera objects and reclaim memory.

    ct.isLoaded(true);

    classList = {'ct.Interface', 'ct.Kinetics', 'ct.Mixture', 'ct.ThermoPhase', ...
                 'ct.Transport', 'ct.Solution', 'ct.Func1', 'ct.Domain1D', ...
                 'ct.Connector', 'ct.ReactorBase', 'ct.ReactorNet'};

    varList = evalin('base', 'whos');

    for i = 1:length(varList)
        for j = 1:length(classList)
            if isa(evalin('base', ['eval("', varList(i).name, '")']), ...
                  classList{j});
                evalin('base', ['clear ', varList(i).name]);
                break
            end
        end
    end
end
