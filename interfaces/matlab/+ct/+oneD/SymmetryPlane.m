classdef SymmetryPlane < ct.oneD.Boundary
    % Create a symmetry plane domain. ::
    %
    %     >> m = ct.oneD.SymmetryPlane(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`ct.oneD.SymmetryPlane`.

    methods

        function obj = SymmetryPlane(phase, name)
            % Constructor
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "symmetry-plane"
            end

            obj@ct.oneD.Boundary('symmetry-plane', phase, name);

        end

    end

end
