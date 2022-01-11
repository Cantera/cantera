function r = IdealGasReactor(contents)
    % Create a ideal gas reactor object.
    % The governing equations are specialized for the ideal gas equation of
    % state (and do not work correctly with other thermodynamic models).
    %
    %: param contents:
    %    Contents of the reactor of class 'Solution'.

    if nargin == 0
        contents = 0;
    end
    r = Reactor(contents, 'IdealGasReactor');
end
