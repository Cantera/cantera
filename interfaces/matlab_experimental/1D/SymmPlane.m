function m = SymmPlane(id)
    % Create a symmetry plane domain.
    % m = SymmPlane(id)
    %
    % :param id:
    %     String ID of the symmetry plane.
    % :return:
    %     Instance of class :mat:func:`Domain1D` representing a symmetry
    %     plane.
    %
    m = Domain1D(Symm1D);
    if nargin == 0
        m.setID('symmetry_plane');
    else
        m.setID(id);
    end
end
