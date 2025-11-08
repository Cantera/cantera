classdef Surface1D < ct.Boundary1D
    % Create a surface domain. ::
    %
    %     >> m = ct.Surface1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of surface

    methods

        %% Surface1D Class Constructor

        function m = Surface1D(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "surface"
            end

            m@ct.Boundary1D('surface', phase, name);

        end

    end

end
