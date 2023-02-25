classdef Surface < Domain1D
    % Create a surface domain. ::
    %
    %     >> m = Surface(id, surface_mech)
    %
    % :param id:
    %     String ID of surface
    % :param surface_mech:
    %     Instance of class :mat:class:`Interface` defining
    %     the surface reaction mechanism to be used. Optional.
    % :return:
    %     Instance of class :mat:class:`Surface` representing a
    %     non-reacting or reacting surface.

    methods

        function m = Surface(id, surface_mech)
            % Constructor

            if nargin < 2
                param = {'Surf1D'};
            else
                param = {'ReactingSurface', surface_mech};
            end

            m@Domain1D(param{:});
            if nargin == 0
                m.setID('surface');
            elseif nargin == 1
                m.setID(id);
            end

        end

    end

end
