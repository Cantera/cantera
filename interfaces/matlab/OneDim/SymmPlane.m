classdef SymmPlane < Boundary1D
    % Create a symmetry plane domain. ::
    %
    %     >> m = SymmPlane(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`SymmPlane`.

    methods

        function m = SymmPlane(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "symmetry-plane"
            end

            m@Boundary1D('symmetry-plane', phase, name);

        end

    end

end
