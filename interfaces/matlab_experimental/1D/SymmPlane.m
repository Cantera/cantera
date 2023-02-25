classdef SymmPlane < Domain1D
    % Create a symmetry plane domain. ::
    %
    %     >> m = SymmPlane(id)
    %
    % :param id:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:class:`SymmPlane`.

    methods

        function m = SymmPlane(id)
            % Constructor

            m@Domain1D('Symm1D');

            if nargin == 0
                m.setID('symmetry_plane');
            else
                m.setID(id);
            end

        end

    end

end
