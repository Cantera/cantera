classdef Surface < Boundary1D
    % Create a surface domain. ::
    %
    %     >> m = Surface(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of surface
    % :return:
    %     Instance of class :mat:class:`Surface` representing a
    %     non-reacting surface.

    methods

        %% Surface Class Constructor

        function m = Surface(phase, name)
            arguments
                phase (1,1) Solution
                name (1,1) string = "surface"
            end

            m@Boundary1D('surface', phase, name);

        end

    end

end
