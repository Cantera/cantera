classdef Surface < ct.oneD.Boundary
    % Create a surface domain. ::
    %
    %     >> m = ct.oneD.Surface(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of surface

    methods

        %% Surface Class Constructor

        function obj = Surface(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "surface"
            end

            obj@ct.oneD.Boundary('surface', phase, name);

        end

    end

end
