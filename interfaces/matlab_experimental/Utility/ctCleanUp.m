function ctCleanUp()
    % Delete all stored Cantera objects and reclaim memory.

    ctIsLoaded;

    classList = {'Interface', 'Kinetics', 'Mixture', 'ThermoPhase', 'Transport', ...
                 'Solution', 'Func1', 'Domain1D', 'Connector', 'ReactorBase', 'ReactorNet'};

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
