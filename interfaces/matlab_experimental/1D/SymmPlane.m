function m = SymmPlane(id)
    % Create an symmetry plane domain.
    m = Domain1D(4);
    if nargin == 0
        m.setID('symmetry_plane');
    else
        m.setID(id);
    end
end
