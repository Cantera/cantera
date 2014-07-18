function n = nAtoms(tp, k, m)
% NATOMS - Number of atoms of element m in species k.

if nargin == 3
    if ischar(m)
        m = elementIndex(tp, m);
    end
    if ischar(k)
        k = speciesIndex(tp, k);
    end
    n = phase_get(tp.tp_id,14, k, m);
else
    error('nAtoms expects three input arguments.')
end
