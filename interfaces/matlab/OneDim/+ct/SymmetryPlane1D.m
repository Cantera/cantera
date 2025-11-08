classdef SymmetryPlane1D < ct.Boundary1D
    % Create a symmetry plane domain. ::
    %
    %     >> m = ct.SymmetryPlane1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`ct.SymmetryPlane1D`.

    methods

        function m = SymmetryPlane1D(phase, name)
            % Constructor
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "symmetry-plane"
            end

            m@ct.Boundary1D('symmetry-plane', phase, name);

        end

    end

end
