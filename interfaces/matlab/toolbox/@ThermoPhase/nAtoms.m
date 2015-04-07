function n = nAtoms(a,k,m)
% NATOMS - Number of atoms of element m in species k.
if nargin == 3
    if ischar(m)
        m = elementIndex(a,m);
    end
    if ischar(k)
        k = speciesIndex(a,k);
    end
    n = phase_get(a.tp_id,14,k,m);
else
    error('usage: nAtoms(phase, k, m)')
end
