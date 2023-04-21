classdef SymmPlane < Boundary1D
    % Create a symmetry plane domain. ::
    %
    %     >> m = SymmPlane(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`SymmPlane`.

    methods

        function m = SymmPlane(phase, id)
            % Constructor

            if nargin < 2
                id = 'symmetry-plane';
            end

            m@Boundary1D('symmetry-plane', phase, id);

        end

    end

end
