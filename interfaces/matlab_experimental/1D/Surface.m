function m = Surface(id, surface_mech)
    % Create a surface domain.
    % m = Surface(id, surface_mech)
    %
    % :param id:
    %     String ID of surface
    % :param surface_mech:
    %     Instance of class :mat:func:`Interface` defining
    %     the surface reaction mechanism to be used. Optional.
    % :return:
    %     Instance of class :mat:func:`Domain1D` representing a
    %     non-reacting or reacting surface.
    %
    if nargin < 2
        m = Domain1D('Surf1D');
        if nargin == 0
            m.setID('surface');
        elseif nargin == 1
            m.setID(id);
        end
    else
        m = Domain1D('ReactingSurface', surface_mech);
        m.setID(id);
    end
end
