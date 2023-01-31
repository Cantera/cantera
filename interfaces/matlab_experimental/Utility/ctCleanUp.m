function ctCleanUp()
    % Delete all stored Cantera objects and reclaim memory.
    %
    ctIsLoaded;

    classList = {'Domain1D', 'Sim1D', 'Func', 'Kinetics', ...
    'Interface', 'Mixture', 'Solution', 'FlowDevice', 'Reactor', ...
    'Wall', 'ReactorNet', 'ReactorSurface', 'ThermoPhase', 'Transport'};

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
