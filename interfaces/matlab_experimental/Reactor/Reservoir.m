function r = Reservoir(contents)
    % Create a reservoir object.
    % An instance of class 'Reactor' configured so that its intensive state
    % is constant in time. A reservoir may be thought of as inifinite in
    % extent, perfectly mixed, and non-reacting, so that fluid may be
    % extracted or added without changing the composition or thermodynamic
    % state. Note that even if the reaction mechanism associated with the
    % fluid in the reactor defines reactions, they are disabled within
    % reservoirs.
    %
    %: param contents:
    %    Contents of the reactor of class 'Solution'.

    if nargin == 0
        contents = 0;
    end
    r = Reactor(contents, 'Reservoir');
end
