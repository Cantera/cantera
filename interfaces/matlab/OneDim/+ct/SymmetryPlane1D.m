classdef SymmetryPlane1D < Boundary1D
    % Create a symmetry plane domain. ::
    %
    %     >> m = SymmetryPlane1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`SymmetryPlane1D`.

    methods

        function m = SymmetryPlane1D(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "symmetry-plane"
            end

            m@Boundary1D('symmetry-plane', phase, name);

        end

    end

end
