classdef Surface < Domain1D
    % Create a surface domain.
    %
    % m = Surface(id, surface_mech)
    %
    % :param id:
    %     String ID of surface
    % :param surface_mech:
    %     Instance of class :mat:class:`Interface` defining
    %     the surface reaction mechanism to be used. Optional.
    % :return:
    %     Instance of class :mat:class:`Surface` representing a
    %     non-reacting or reacting surface.
    %

    methods

        % Constructor
        function m = Surface(id, surface_mech)

            if nargin > 1
                m = m@Domain1D('ReactingSurface', surface_mech);
                m.setID(id);
                return
            end

            m = m@Domain1D('Surf1D');

            if nargin == 0
                m.setID('surface');
            elseif nargin == 1
                m.setID(id);
            end

        end

    end

end
