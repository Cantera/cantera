function r = FlowReactor(contents)
    % Create a flow reactor object.
    % A reactor representing adiabatic plug flow in a constant-area duct.
    %
    %: param contents:
    %    Contents of the reactor of class 'Solution'.

    if nargin == 0
        contents = 0;
    end
    r = Reactor(contents, 'FlowReactor');
end
