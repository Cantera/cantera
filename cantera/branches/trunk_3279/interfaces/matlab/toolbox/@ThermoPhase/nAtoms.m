function n = nAtoms(tp,k,m)
% NATOMS  Get the number of atoms of an element in a species.
% n = nAtoms(tp,k,m)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param k:
%     String species name or integer species number
% :param m:
%     String element name or integer element number
% :return:
%     Number of atoms of element ``m`` in species ``k``.
%

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
