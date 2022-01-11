function m = Surface(id, surface_mech)
    % Create an surface domain.
    % :param surface_mech
    %    Instance of class 'Interface' defining the surface reaction
    %    mechanism to be used. Optional.
    if nargin < 2
        m = Domain1D(3);
        if nargin == 0
            m.setID('surface');
        elseif nargin == 1
            m.setID(id);
        end
    else
        m = Domain1D(6, surface_mech);
        m.setID(id);
    end
end
