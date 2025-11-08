classdef Surface1D < Boundary1D
    % Create a surface domain. ::
    %
    %     >> m = Surface1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of surface

    methods

        %% Surface1D Class Constructor

        function m = Surface1D(phase, name)
            arguments
                phase (1,1) Solution
                name (1,1) string = "surface"
            end

            m@Boundary1D('surface', phase, name);

        end

    end

end
