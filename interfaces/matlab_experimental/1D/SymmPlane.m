function m = SymmPlane(id)
    % Create an symmetry plane domain.
    m = Domain1D(Symm1D);
    if nargin == 0
        m.setID('symmetry_plane');
    else
        m.setID(id);
    end
end
